// DIFFCP SOURCE: diffcp/cpp/src/lsqr.cpp (Eigen LSQR wrapper).
//
// Phase 2 scaffolding. The full LSQR-on-LinearOperator implementation
// lands in a follow-on commit; this file establishes the
// RcppEigen LinkingTo plumbing so the package builds cleanly with C++
// from day one.

#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Rcpp;

// [[Rcpp::export]]
double diffcp_dot(Eigen::Map<Eigen::VectorXd> x,
                  Eigen::Map<Eigen::VectorXd> y) {
  // Smoke-test entry point: confirms the C++ build links against
  // RcppEigen and that vectors round-trip from R.
  return x.dot(y);
}
