## PSD-cone tests, ported from
##   /Users/naras/research/cvxr/new_design/diffing/diffcp/tests/test_clarabel_psd.py
##
## Each Python test builds a problem in cvxpy, canonicalises to SCS form
## via `scs_data_from_cvxpy_problem`, then solves via `diffcp.solve_and_derivative`
## with both SCS and Clarabel and asserts:
##   * obj_solver ~= cvxpy_obj at atol=1e-4
##   * x_scs ~= x_clarabel at atol=1e-4
##   * (some tests) s_scs ~= s_clarabel at atol=1e-4
##   * test_constraint_satisfaction: A @ x + s ~= b at atol=1e-5
##   * test_derivative_lsqr_mode: DT() runs in lsqr mode and is finite
##
## We don't have cvxpy in R; the canonicalised (A, b, c, cone_dict)
## and `cvxpy_obj` reference come from inst/python_parity/psd.py and
## are pinned here as numeric literals.

skip_if_not_installed("clarabel")
skip_if_not_installed("scs")

## Helper: build sparse A from CSC triplets emitted by the parity script
## (1-based row indices and 1-based indptr).
.csc <- function(A_i, A_p, A_x, dims) {
  Matrix::sparseMatrix(i = A_i, p = A_p - 1L, x = A_x, dims = dims, index1 = TRUE)
}

## Helper: solve via both Clarabel and SCS, return list of solutions.
.solve_both <- function(A, b, c, cone_dict, ctrl_clarabel, ctrl_scs) {
  res_cla <- do.call(diffcp::solve_only,
                     c(list(A, b, c, cone_dict, solve_method = "Clarabel"),
                       ctrl_clarabel))
  res_scs <- do.call(diffcp::solve_only,
                     c(list(A, b, c, cone_dict, solve_method = "SCS"),
                       ctrl_scs))
  list(clarabel = res_cla, scs = res_scs)
}

# ---- Tighter tolerances than Python's atol=1e-4 ------------------
ctrl_cla <- list(tol_gap_abs = 1e-9, tol_gap_rel = 1e-9, tol_feas = 1e-9)
ctrl_scs <- list(eps_abs = 1e-9, eps_rel = 1e-9, max_iters = 100000L)

# ---- Port of test_multiple_psd_cones_objective_match -------------
test_that("multiple PSD cones: objective and x match between SCS and Clarabel", {
  A_i <- c(1, 10, 1, 11, 1, 12, 1, 13, 1, 14, 1, 15, 2, 4, 5, 6, 7, 8, 9,
           3, 4, 5, 6, 7, 8, 9)
  A_p <- c(1, 3, 5, 7, 9, 11, 13, 20, 27)
  A_x <- c(1, -1, 4, -1.414213562373095, 6, -1.414213562373095,
           4, -1, 10, -1.414213562373095, 6, -1, -1, -1,
           -2.82842712474619, -4.242640687119286, -4, -7.071067811865476, -6,
           -1, -7, -11.31370849898476, -12.72792206135786, -10,
           -15.55634918610405, -12)
  A <- .csc(A_i, A_p, A_x, c(15, 8))
  b <- c(1, rep(0, 14))
  c <- c(1, 0, 0, 1, 0, 1, 1, 1)
  cone_dict <- list(l = 2L, s = c(3L, 3L), z = 1L)
  obj_ref <- 0.08814600004013824

  sol <- .solve_both(A, b, c, cone_dict, ctrl_cla, ctrl_scs)

  obj_cla <- sum(c * sol$clarabel$x)
  obj_scs <- sum(c * sol$scs$x)
  expect_lt(abs(obj_cla - obj_ref), 1e-4)
  expect_lt(abs(obj_scs - obj_ref), 1e-4)
  expect_lt(max(abs(sol$clarabel$x - sol$scs$x)), 1e-4)
  expect_lt(max(abs(sol$clarabel$s - sol$scs$s)), 1e-4)
})

# ---- Port of test_single_psd_cone --------------------------------
test_that("single PSD cone: objective and x match between SCS and Clarabel", {
  A_i <- c(1, 2, 3, 4, 1, 5, 6, 1, 7)
  A_p <- c(1, 3, 4, 5, 7, 8, 10)
  A_x <- c(1, -1, -1.414213562373095, -1.414213562373095,
           1, -1, -1.414213562373095, 1, -1)
  A <- .csc(A_i, A_p, A_x, c(7, 6))
  b <- c(1, rep(0, 6))
  c <- c(1, 0, 0, 1, 0, 1)
  cone_dict <- list(s = c(3L), z = 1L)
  obj_ref <- 1.000000000000004

  sol <- .solve_both(A, b, c, cone_dict, ctrl_cla, ctrl_scs)
  obj_cla <- sum(c * sol$clarabel$x)
  obj_scs <- sum(c * sol$scs$x)
  expect_lt(abs(obj_cla - obj_ref), 1e-4)
  expect_lt(abs(obj_scs - obj_ref), 1e-4)
  expect_lt(max(abs(sol$clarabel$x - sol$scs$x)), 1e-4)
  expect_lt(max(abs(sol$clarabel$s - sol$scs$s)), 1e-4)
})

# ---- Port of test_mixed_cones ------------------------------------
test_that("mixed cones (z + l + PSD): objective and x match between SCS and Clarabel", {
  A_i <- c(1, 4, 1, 5, 1, 6, 2, 3)
  A_p <- c(1, 3, 5, 7, 9)
  A_x <- c(1, -1, 1, -1.414213562373095, 2, -1, -1, 1)
  A <- .csc(A_i, A_p, A_x, c(6, 4))
  b <- c(1, 0, 5, 0, 0, 0)
  c <- c(1, 0, 1, 1)
  cone_dict <- list(l = 2L, s = c(2L), z = 1L)
  obj_ref <- 0.4530818398464256

  sol <- .solve_both(A, b, c, cone_dict, ctrl_cla, ctrl_scs)
  obj_cla <- sum(c * sol$clarabel$x)
  obj_scs <- sum(c * sol$scs$x)
  expect_lt(abs(obj_cla - obj_ref), 1e-4)
  expect_lt(abs(obj_scs - obj_ref), 1e-4)
  expect_lt(max(abs(sol$clarabel$x - sol$scs$x)), 1e-4)
})

# ---- Port of test_constraint_satisfaction ------------------------
test_that("Clarabel solution satisfies A @ x + s = b for multi-PSD problem", {
  A_i <- c(1, 10, 1, 11, 1, 12, 1, 13, 1, 14, 1, 15, 2, 4, 5, 6, 7, 8, 9,
           3, 4, 5, 6, 7, 8, 9)
  A_p <- c(1, 3, 5, 7, 9, 11, 13, 20, 27)
  A_x <- c(1, -1, 4, -1.414213562373095, 6, -1.414213562373095,
           4, -1, 10, -1.414213562373095, 6, -1, -1, -1,
           -2.82842712474619, -4.242640687119286, -4, -7.071067811865476, -6,
           -1, -7, -11.31370849898476, -12.72792206135786, -10,
           -15.55634918610405, -12)
  A <- .csc(A_i, A_p, A_x, c(15, 8))
  b <- c(1, rep(0, 14))
  c <- c(1, 0, 0, 1, 0, 1, 1, 1)
  cone_dict <- list(l = 2L, s = c(3L, 3L), z = 1L)

  res <- do.call(diffcp::solve_only,
                 c(list(A, b, c, cone_dict, solve_method = "Clarabel"),
                   ctrl_cla))
  residual <- as.numeric(A %*% res$x + res$s - b)
  expect_lt(max(abs(residual)), 1e-5)
})

# ---- Port of test_derivative_lsqr_mode ---------------------------
test_that("Clarabel + lsqr mode: DT runs and returns finite values (PSD)", {
  A_i <- c(1, 2, 3, 1, 4)
  A_p <- c(1, 3, 4, 6)
  A_x <- c(1, -1, -1.414213562373095, 1, -1)
  A <- .csc(A_i, A_p, A_x, c(4, 3))
  b <- c(1, 0, 0, 0)
  c <- c(1, 0.6, 2)
  cone_dict <- list(s = c(2L), z = 1L)

  res <- do.call(diffcp::solve_and_derivative,
                 c(list(A, b, c, cone_dict,
                        solve_method = "Clarabel", mode = "lsqr"),
                   ctrl_cla))

  set.seed(42)
  dx <- rnorm(length(res$x)) * 0.01
  dy <- rnorm(length(res$y)) * 0.01
  ds <- rnorm(length(res$s)) * 0.01

  da <- res$DT(dx, dy, ds)
  expect_true(all(is.finite(da$dc)), info = "dc has non-finite values")
  expect_true(all(is.finite(da$db)), info = "db has non-finite values")
  expect_true(all(is.finite(methods::as(da$dA, "CsparseMatrix")@x)),
              info = "dA has non-finite values")
})

# ---- Port of test_psd_permutation_logic --------------------------
test_that("PSD permutation logic: random C agrees with cvxpy reference", {
  A_i <- c(1, 2, 3, 4, 1, 5, 6, 1, 7)
  A_p <- c(1, 3, 4, 5, 7, 8, 10)
  A_x <- c(1, -1, -1.414213562373095, -1.414213562373095,
           1, -1, -1.414213562373095, 1, -1)
  A <- .csc(A_i, A_p, A_x, c(7, 6))
  b <- c(1, rep(0, 6))
  c <- c(4.229934512253635, 7.487741429710518, 3.028829625895845,
         9.464447231296065, 3.894503620788191, 0.9362311369855318)
  cone_dict <- list(s = c(3L), z = 1L)
  obj_ref <- 0.306530700912452

  res <- do.call(diffcp::solve_only,
                 c(list(A, b, c, cone_dict, solve_method = "Clarabel"),
                   ctrl_cla))
  obj_cla <- sum(c * res$x)
  expect_lt(abs(obj_cla - obj_ref), 1e-4)
})
