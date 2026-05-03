## Forward QP solve via Clarabel and SCS, mirroring how cvxpylayers
## calls `solve_only(P=P)` in Python diffcp.
##
## Python diffcp's QP support is forward-only: solve_and_derivative
## with P != NULL and mode in {"dense", "lsqr", "lsmr"} explicitly
## errors and points the user to lpgd mode (cone_program.py:637-639).
## We mirror both behaviours here:
##   * solve_only(P) returns x/y/s matching Python's reference for
##     the same problem.
##   * solve_and_derivative(P) errors with Python's exact message
##     pointing to lpgd.
##
## Reference numerics from inst/python_parity/qp.py.

skip_if_not_installed("clarabel")
skip_if_not_installed("scs")

.csc <- function(A_i, A_p, A_x, dims) {
  Matrix::sparseMatrix(i = A_i, p = A_p - 1L, x = A_x, dims = dims, index1 = TRUE)
}

A <- .csc(c(1, 2, 1, 3, 1, 4), c(1, 3, 5, 7),
          c(1, -1, 1, -1, 1, -1), c(4, 3))
P <- .csc(c(1, 2, 1, 2, 3, 2, 3), c(1, 3, 6, 8),
          c(2, 0.5, 0.5, 2, 0.5, 0.5, 2), c(3, 3))
b <- c(1, 0, 0, 0)
c <- c(-1, 0.5, 1)
cone_dict <- list(z = 1L, l = 3L)

test_that("solve_only forward QP matches Python diffcp (Clarabel)", {
  res <- diffcp::solve_only(A, b, c, cone_dict, P = P,
                            solve_method = "Clarabel",
                            tol_gap_abs = 1e-12,
                            tol_gap_rel = 1e-12,
                            tol_feas    = 1e-12)
  expect_equal(res$x,
               c(0.9999997064006334, 1.590246952183832e-07,
                 1.345746671188647e-07),
               tolerance = 1e-7)
  expect_equal(res$y,
               c(-0.9999994923135322, 9.167811022181554e-14,
                  7.4622350325403e-07,  8.56348145174723e-07),
               tolerance = 1e-7)
  expect_equal(res$s,
               c(0, 0.9999997064006332, 1.590246889842621e-07,
                  1.345746599683592e-07),
               tolerance = 1e-7)
})

test_that("solve_only forward QP matches Python diffcp (SCS)", {
  res <- diffcp::solve_only(A, b, c, cone_dict, P = P,
                            solve_method = "SCS",
                            eps_abs = 1e-12, eps_rel = 1e-12,
                            max_iters = 100000L)
  expect_equal(res$x,
               c(1, 5.051524543201496e-17, -2.971459647020941e-17),
               tolerance = 1e-7)
  expect_equal(res$s,
               c(0, 1, 0, 0),
               tolerance = 1e-7)
})

test_that("solve_and_derivative aborts on QP with the upstream message", {
  ## Python diffcp/cone_program.py:637-639 raises:
  ##   "Dense, lsqr, and lsmr modes currently do not support quadratic
  ##    objectives. Consider switching to 'lpgd' mode."
  ## We mirror that exact text (lpgd is not yet ported in our R port).
  expect_error(
    diffcp::solve_and_derivative(A, b, c, cone_dict, P = P,
                                 mode = "dense"),
    "do not support quadratic objectives"
  )
  expect_error(
    diffcp::solve_and_derivative(A, b, c, cone_dict, P = P,
                                 mode = "lsqr"),
    "do not support quadratic objectives"
  )
})
