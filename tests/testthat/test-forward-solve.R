## Forward-solve parity vs Python `diffcp.cone_program.solve_only`
## (Clarabel backend). Expected (x, y, s) values were generated via:
##
##     uv run python inst/python_parity/forward_solve.py
##
## and pinned here as numeric literals so the test does not depend on
## Python at run time.

skip_if_not_installed("clarabel")

# -- LP with one equality + nonneg orthant -------------------------
test_that("solve_only matches Python diffcp on a small LP", {
  ## minimise c^T x  s.t. 1^T x = 1, x >= 0   (n = 3)
  A_data <- c(1, -1, 0, 0,
              1,  0, -1, 0,
              1,  0, 0, -1)
  A <- Matrix::Matrix(matrix(A_data, nrow = 4, ncol = 3), sparse = TRUE)
  b <- c(1, 0, 0, 0)
  c <- c(1, 2, 3)
  cone_dict <- list(z = 1L, l = 3L)

  res <- diffcp::solve_only(A, b, c, cone_dict,
                            tol_gap_abs = 1e-10,
                            tol_gap_rel = 1e-10,
                            tol_feas    = 1e-10)

  expect_equal(res$x,
               c(0.999999999979, 1.56223305592e-11, 5.32985362563e-12),
               tolerance = 1e-7)
  expect_equal(res$y,
               c(-1.00000000003, 6.99807870964e-12,
                 1.00000000001, 2.00000000001),
               tolerance = 1e-7)
  expect_equal(res$s,
               c(0, 0.99999999999, 2.64458841362e-11, 1.61534072027e-11),
               tolerance = 1e-7)
})

# -- SOC: minimise t s.t. ||x||_2 <= t, sum(x)=1 -------------------
test_that("solve_only matches Python diffcp on a small SOC program", {
  A_data <- c(0, -1, 0, 0, 0,
              1,  0, -1, 0, 0,
              1,  0, 0, -1, 0,
              1,  0, 0, 0, -1)
  A <- Matrix::Matrix(matrix(A_data, nrow = 5, ncol = 4), sparse = TRUE)
  b <- c(1, 0, 0, 0, 0)
  c <- c(1, 0, 0, 0)
  cone_dict <- list(z = 1L, q = 4L)

  res <- diffcp::solve_only(A, b, c, cone_dict,
                            tol_gap_abs = 1e-10,
                            tol_gap_rel = 1e-10,
                            tol_feas    = 1e-10)

  expect_equal(res$x,
               c(0.577350269148, 0.333333333333,
                 0.333333333333, 0.333333333333),
               tolerance = 1e-7)
  expect_equal(res$y,
               c(-0.577350269158, 1,
                 -0.577350269158, -0.577350269158, -0.577350269158),
               tolerance = 1e-7)
  expect_equal(res$s,
               c(0, 0.577350269229,
                 0.333333333333, 0.333333333333, 0.333333333333),
               tolerance = 1e-7)
})

# -- Least-1-norm with simplex constraints (m=4, n=2) --------------
test_that("solve_only matches Python diffcp on a least-1-norm fit", {
  A_data <- c(
    1, -1, 0,
    -1.76405234597, -0.978737984106, -1.86755799015, -0.950088417526,
     1.76405234597,  0.978737984106,  1.86755799015,  0.950088417526,
    1,  0, -1,
    -0.400157208367, -2.2408931992, 0.977277879876, 0.151357208298,
     0.400157208367,  2.2408931992, -0.977277879876, -0.151357208298,
    0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0,
    0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0,
    0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0,
    0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1
  )
  A <- Matrix::Matrix(matrix(A_data, nrow = 11, ncol = 6), sparse = TRUE)
  b <- c(1, 0, 0,
         0.103218851794, -0.410598501938, -0.144043571161, -1.45427350696,
        -0.103218851794,  0.410598501938,  0.144043571161,  1.45427350696)
  c <- c(0, 0, 1, 1, 1, 1)
  cone_dict <- list(z = 1L, l = 10L)

  res <- diffcp::solve_only(A, b, c, cone_dict,
                            tol_gap_abs = 1e-10,
                            tol_gap_rel = 1e-10,
                            tol_feas    = 1e-10)

  expect_equal(res$x,
               c(0.394160332091, 0.605839667909,
                 1.04096942035, 1.33280317826,
                 2.71078680966e-10, 1.1714845415),
               tolerance = 1e-6)
  expect_equal(res$y,
               c(-2.44898173259,
                 9.86831015808e-10, 5.27248031171e-11, 9.65403455408e-11,
                 1.00341548699e-10, 0.324294660638,
                 1.00000000015, 1.00000000006, 1.00000000003,
                 0.675705339506, 4.35026702385e-11),
               tolerance = 1e-6)
  expect_equal(res$s,
               c(0,
                 0.394160332234, 0.60583966804,
                 2.08193884114, 2.66560635708,
                 7.75671769009e-10, 8.63971826889e-11,
                 6.73624518845e-11, 2.38683178674e-11,
                 3.01504984204e-10, 2.3429690833),
               tolerance = 1e-6)
})

test_that("solve_only rejects warm_start for Clarabel", {
  A <- Matrix::Matrix(matrix(c(1, -1, 0, 0, 1, 0, -1, 0, 1, 0, 0, -1),
                             nrow = 4, ncol = 3), sparse = TRUE)
  b <- c(1, 0, 0, 0)
  c <- c(1, 2, 3)
  expect_error(
    diffcp::solve_only(A, b, c, list(z = 1L, l = 3L),
                       warm_start = list(c(0, 0, 0), c(0, 0, 0, 0), c(0, 0, 0, 0))),
    "warm.start"
  )
})

test_that("solve_and_derivative still errors with phase-pointer", {
  A <- Matrix::Matrix(matrix(c(1, -1, 0, 0, 1, 0, -1, 0, 1, 0, 0, -1),
                             nrow = 4, ncol = 3), sparse = TRUE)
  b <- c(1, 0, 0, 0)
  c <- c(1, 2, 3)
  expect_error(
    diffcp::solve_and_derivative(A, b, c, list(z = 1L, l = 3L)),
    "not yet implemented"
  )
})
