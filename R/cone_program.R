## CVXPY SOURCE: diffcp/cone_program.py
##
## The main solve_and_derivative entry point. Phase 2 scaffolding;
## the full implementation lands incrementally as cones (Phase 2a),
## the M operator (Phase 2b), and LSQR (Phase 2c) come online.

#' Solve a cone program and return forward / adjoint derivative operators
#'
#' Mirrors `diffcp.cone_program.solve_and_derivative`. Solves
#'   minimise  c^T x        s.t.  A x + s = b,  s in K
#' (with optional QP `0.5 x^T P x` term) and returns the optimal
#' `(x, y, s)` together with closures `D` (forward) and `DT` (adjoint)
#' that map perturbations of `(A, b, c, [P])` to perturbations of
#' `(x, y, s)` and vice versa.
#'
#' @param A A sparse `dgCMatrix` constraint matrix.
#' @param b A numeric offset vector.
#' @param c A numeric objective coefficient vector.
#' @param cones A cone list as produced by `parse_cone_dict`.
#' @param P Optional sparse `dgCMatrix` for QP objective; default `NULL`.
#' @param solver Solver name to call for the forward solve (e.g.
#'   `"Clarabel"`, `"SCS"`); default `"Clarabel"`.
#' @param mode Differentiation mode: `"lsqr"` (default), `"dense"`,
#'   or `"lpgd"`.
#' @param warm_start Optional warm-start `(x, y, s)` triple.
#' @param ... Additional solver kwargs.
#' @returns A list with elements `x`, `y`, `s`, `D`, `DT`.
#' @export
solve_and_derivative <- function(A, b, c, cones,
                                 P = NULL,
                                 solver = "Clarabel",
                                 mode = "lsqr",
                                 warm_start = NULL,
                                 ...) {
  cli::cli_abort(c(
    "{.fn solve_and_derivative} not yet implemented (Phase 2 scaffolding).",
    "i" = "Cone projections (Zero, Nonneg, SOC) and their Jacobians are in place.",
    "i" = "The forward solve, M operator, and LSQR-based derivative arrive in subsequent commits per {.file notes/diffcp_r_plan.md}."
  ))
}
