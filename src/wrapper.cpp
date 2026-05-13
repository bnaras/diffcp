// DIFFCP SOURCE: diffcp/cpp/src/wrapper.cpp @ v1.1.8 (6a78143) (replaced pybind11 with Rcpp).
// Upstream pin: see inst/UPSTREAM.dcf / upstream_info().
//
// Rcpp-callable entry points to the C++ derivative-machinery. The R
// side hands us:
//   * cones encoded as a List of List(name = "z"/"l"/..., sizes =
//     IntegerVector). We translate to std::vector<Cone> on entry.
//   * sparse matrices as Eigen::SparseMatrix<double> via RcppEigen.
//   * dense matrices/vectors as Eigen::MatrixXd / Eigen::VectorXd.
//
// All the heavy work (LinearOperator composition, Eigen::LDLT solves,
// matrix-free LSQR) happens entirely in C++; the R side just builds
// inputs and consumes outputs.

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#include "cones.h"
#include "deriv.h"
#include "linop.h"
#include "lsqr.h"

using Rcpp::as;
using Rcpp::List;
using Rcpp::IntegerVector;

// -- Cone-dict translation ------------------------------------------
//
// R-side schema: List, where each element is List(name = chr,
// sizes = int).  The name keys mirror SCS / Clarabel ("z", "l", "q",
// "s", "ep", "ed").
static ConeType cone_type_from_name(const std::string &name) {
  if (name == "z" || name == "f") return ZERO;
  if (name == "l")                return POS;
  if (name == "q")                return SOC;
  if (name == "s")                return PSD;
  if (name == "ep")               return EXP;
  if (name == "ed")               return EXP_DUAL;
  Rcpp::stop("Unknown cone name: '%s'", name);
}

static std::vector<Cone> parse_cones(const List &cones_R) {
  std::vector<Cone> out;
  out.reserve(cones_R.size());
  for (int i = 0; i < cones_R.size(); ++i) {
    List entry = cones_R[i];
    std::string name = as<std::string>(entry["name"]);
    IntegerVector sizes_R = as<IntegerVector>(entry["sizes"]);
    std::vector<int> sizes(sizes_R.begin(), sizes_R.end());
    out.emplace_back(cone_type_from_name(name), sizes);
  }
  return out;
}

// -- Cone projection Jacobian, dense --------------------------------

// [[Rcpp::export]]
Eigen::MatrixXd cpp_dprojection_dense(const Eigen::VectorXd &x,
                                      const List &cones_R, bool dual) {
  return dprojection_dense(x, parse_cones(cones_R), dual);
}

// -- Dense M = (Q - I) D pi(z) + I ----------------------------------
//
// Q must be passed dense.  This matches diffcp/cpp/src/deriv.cpp's
// M_dense signature; the Python source materialises Q with
// `Q.todense()` before the call.

// [[Rcpp::export]]
Eigen::MatrixXd cpp_M_dense(const Eigen::MatrixXd &Q_dense,
                            const List &cones_R,
                            const Eigen::VectorXd &u,
                            const Eigen::VectorXd &v,
                            double w) {
  return M_dense(Q_dense, parse_cones(cones_R), u, v, w);
}

// -- Dense derivative / adjoint solves via Eigen LDLT ---------------
//
// Both functions form the appropriate Gram matrix on the C++ side and
// factor it with LDLT, exactly as in deriv.cpp:53-63.

// [[Rcpp::export]]
Eigen::VectorXd cpp_solve_derivative_dense(const Eigen::MatrixXd &M,
                                           const Eigen::MatrixXd &MT,
                                           const Eigen::VectorXd &rhs) {
  return _solve_derivative_dense(M, MT, rhs);
}

// [[Rcpp::export]]
Eigen::VectorXd cpp_solve_adjoint_derivative_dense(const Eigen::MatrixXd &M,
                                                   const Eigen::MatrixXd &MT,
                                                   const Eigen::VectorXd &dz) {
  return _solve_adjoint_derivative_dense(M, MT, dz);
}

// -- Matrix-free LSQR against the M operator ------------------------
//
// Build the matrix-free M from (Q, cones, u, v, w), then run LSQR
// against either M or M^T depending on `transpose`.

// [[Rcpp::export]]
Rcpp::List cpp_lsqr_M(const Eigen::Map<Eigen::SparseMatrix<double>> &Q,
                      const List &cones_R,
                      const Eigen::VectorXd &u,
                      const Eigen::VectorXd &v,
                      double w,
                      const Eigen::VectorXd &rhs,
                      bool transpose,
                      double damp = 0.0,
                      double atol = 1e-8,
                      double btol = 1e-8,
                      double conlim = 1e8,
                      int iter_lim = -1) {
  std::vector<Cone> cones = parse_cones(cones_R);
  // M_operator returns a LinearOperator wrapping (Q - I) D pi + I.
  // It captures Q by value; we hand it a fresh SparseMatrix so the
  // closure stays alive after this function returns.
  Eigen::SparseMatrix<double> Q_owned = Q;  // copy out of Eigen::Map
  LinearOperator M = M_operator(Q_owned, cones, u, v, w);
  LinearOperator A = transpose ? M.transpose() : M;
  LsqrResult res = lsqr(A, rhs, damp, atol, btol, conlim, iter_lim);
  return Rcpp::List::create(
    Rcpp::Named("solution") = res.x,
    Rcpp::Named("istop")    = res.istop,
    Rcpp::Named("itn")      = res.itn,
    Rcpp::Named("r1norm")   = res.r1norm,
    Rcpp::Named("r2norm")   = res.r2norm,
    Rcpp::Named("anorm")    = res.anorm,
    Rcpp::Named("acond")    = res.acond,
    Rcpp::Named("arnorm")   = res.arnorm,
    Rcpp::Named("xnorm")    = res.xnorm
  );
}

// -- Sparse LSQR (for testing the LSQR routine in isolation) --------

// [[Rcpp::export]]
Rcpp::List cpp_lsqr_sparse(const Eigen::Map<Eigen::SparseMatrix<double>> &A,
                           const Eigen::VectorXd &b,
                           double damp = 0.0,
                           double atol = 1e-8,
                           double btol = 1e-8,
                           double conlim = 1e8,
                           int iter_lim = -1) {
  Eigen::SparseMatrix<double> A_owned = A;
  LsqrResult res = lsqr_sparse(A_owned, b, damp, atol, btol, conlim, iter_lim);
  return Rcpp::List::create(
    Rcpp::Named("solution") = res.x,
    Rcpp::Named("istop")    = res.istop,
    Rcpp::Named("itn")      = res.itn,
    Rcpp::Named("r1norm")   = res.r1norm,
    Rcpp::Named("r2norm")   = res.r2norm,
    Rcpp::Named("anorm")    = res.anorm,
    Rcpp::Named("acond")    = res.acond,
    Rcpp::Named("arnorm")   = res.arnorm,
    Rcpp::Named("xnorm")    = res.xnorm
  );
}

// -- Exponential-cone projection (single 3-vector) ------------------

// [[Rcpp::export]]
Eigen::VectorXd cpp_project_exp_cone(const Eigen::VectorXd &x) {
  if (x.size() != 3) Rcpp::stop("cpp_project_exp_cone requires length-3 input");
  Eigen::Vector3d xv(x[0], x[1], x[2]);
  Eigen::Vector3d out = project_exp_cone(xv);
  Eigen::VectorXd r(3);
  r << out[0], out[1], out[2];
  return r;
}

// [[Rcpp::export]]
bool cpp_in_exp(const Eigen::VectorXd &x) {
  if (x.size() != 3) Rcpp::stop("cpp_in_exp requires length-3 input");
  Eigen::Vector3d xv(x[0], x[1], x[2]);
  return in_exp(xv);
}

// [[Rcpp::export]]
bool cpp_in_exp_dual(const Eigen::VectorXd &x) {
  if (x.size() != 3) Rcpp::stop("cpp_in_exp_dual requires length-3 input");
  Eigen::Vector3d xv(x[0], x[1], x[2]);
  return in_exp_dual(xv);
}
