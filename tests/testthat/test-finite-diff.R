## Finite-difference correctness tests, mirroring Python diffcp's
## `test_solve_and_derivative` (test_clarabel.py:9-43, test_scs.py:9-43).
##
## Idea: solve at (A, b, c), then re-solve at (A + dA, b + db, c + dc)
## for small random (dA, db, dc).  The forward derivative D should
## satisfy  D(dA, db, dc) ~= (x_pert - x, y_pert - y, s_pert - s)
## up to atol=1e-8 (the tolerance Python uses).
##
## For the adjoint: applying DT(c, 0, 0) yields directions
## (dA, db, dc) such that perturbing the problem by 1e-6 * (dA, db, dc)
## changes the objective c^T x by 1e-6 * <(dA, db, dc), (dA, db, dc)>
## (a self-consistency check on the adjoint, also at atol=1e-8).

skip_if_not_installed("clarabel")

test_that("D matches finite-difference (LP, mode='dense')", {
  set.seed(0)
  A_data <- c(1, -1, 0, 0,
              1,  0, -1, 0,
              1,  0, 0, -1)
  A <- Matrix::Matrix(matrix(A_data, nrow = 4, ncol = 3), sparse = TRUE)
  b <- c(1, 0, 0, 0)
  c <- c(1, 2, 3)
  cone_dict <- list(z = 1L, l = 3L)
  ctrl <- list(tol_gap_abs = 1e-12, tol_gap_rel = 1e-12, tol_feas = 1e-12)

  res <- do.call(diffcp::solve_and_derivative,
                 c(list(A, b, c, cone_dict, mode = "dense"), ctrl))

  ## Random perturbation with the same nonzero pattern as A.
  A_csc <- methods::as(A, "CsparseMatrix")
  nz <- Matrix::summary(A_csc)
  dA <- Matrix::sparseMatrix(i = nz$i, j = nz$j,
                             x = rnorm(nrow(nz)) * 1e-6,
                             dims = c(nrow(A), ncol(A)))
  db <- rnorm(length(b)) * 1e-6
  dc <- rnorm(length(c)) * 1e-6

  d <- res$D(dA, db, dc)

  res_pert <- do.call(diffcp::solve_only,
                      c(list(A + dA, b + db, c + dc, cone_dict), ctrl))

  ## Same atol=1e-8 as Python (test_clarabel.py:29-31).  testthat's
  ## tolerance arg is relative; use expect_lt on the max abs residual
  ## to match numpy's `np.testing.assert_allclose(..., atol=1e-8)`.
  expect_lt(max(abs(d$dx - (res_pert$x - res$x))), 1e-8)
  expect_lt(max(abs(d$dy - (res_pert$y - res$y))), 1e-8)
  expect_lt(max(abs(d$ds - (res_pert$s - res$s))), 1e-8)
})

test_that("DT objective sensitivity matches finite-difference (LP, mode='dense')", {
  A_data <- c(1, -1, 0, 0,
              1,  0, -1, 0,
              1,  0, 0, -1)
  A <- Matrix::Matrix(matrix(A_data, nrow = 4, ncol = 3), sparse = TRUE)
  b <- c(1, 0, 0, 0)
  c <- c(1, 2, 3)
  cone_dict <- list(z = 1L, l = 3L)
  ctrl <- list(tol_gap_abs = 1e-12, tol_gap_rel = 1e-12, tol_feas = 1e-12)

  res <- do.call(diffcp::solve_and_derivative,
                 c(list(A, b, c, cone_dict, mode = "dense"), ctrl))
  objective <- sum(c * res$x)

  ## DT(c, 0, 0) gives steepest descent of the objective in (A, b, c).
  da <- res$DT(c, numeric(length(b)), numeric(length(b)))

  eps <- 1e-6
  res_pert <- do.call(diffcp::solve_only,
                      c(list(A + eps * da$dA, b + eps * da$db, c + eps * da$dc,
                             cone_dict), ctrl))
  objective_pert <- sum(c * res_pert$x)

  ## Mirrors test_clarabel.py:41-43:
  ##   assert_allclose(obj_pert - obj,
  ##                   eps * (sum(dA*dA) + db@db + dc@dc),  atol=1e-8)
  dA_sumsq <- sum(da$dA * da$dA)
  expected <- eps * (dA_sumsq + sum(da$db * da$db) + sum(da$dc * da$dc))
  expect_lt(abs(objective_pert - objective - expected), 1e-8)
})

test_that("D matches finite-difference (SOC, mode='dense')", {
  set.seed(1)
  A_data <- c(0, -1, 0, 0, 0,
              1,  0, -1, 0, 0,
              1,  0, 0, -1, 0,
              1,  0, 0, 0, -1)
  A <- Matrix::Matrix(matrix(A_data, nrow = 5, ncol = 4), sparse = TRUE)
  b <- c(1, 0, 0, 0, 0)
  c <- c(1, 0, 0, 0)
  cone_dict <- list(z = 1L, q = 4L)
  ctrl <- list(tol_gap_abs = 1e-12, tol_gap_rel = 1e-12, tol_feas = 1e-12)

  res <- do.call(diffcp::solve_and_derivative,
                 c(list(A, b, c, cone_dict, mode = "dense"), ctrl))

  A_csc <- methods::as(A, "CsparseMatrix")
  nz <- Matrix::summary(A_csc)
  dA <- Matrix::sparseMatrix(i = nz$i, j = nz$j,
                             x = rnorm(nrow(nz)) * 1e-6,
                             dims = c(nrow(A), ncol(A)))
  db <- rnorm(length(b)) * 1e-6
  dc <- rnorm(length(c)) * 1e-6

  d <- res$D(dA, db, dc)
  res_pert <- do.call(diffcp::solve_only,
                      c(list(A + dA, b + db, c + dc, cone_dict), ctrl))

  expect_lt(max(abs(d$dx - (res_pert$x - res$x))), 1e-8)
  expect_lt(max(abs(d$dy - (res_pert$y - res$y))), 1e-8)
  expect_lt(max(abs(d$ds - (res_pert$s - res$s))), 1e-8)
})

test_that("D matches finite-difference (LP, mode='lsqr')", {
  set.seed(2)
  A_data <- c(1, -1, 0, 0,
              1,  0, -1, 0,
              1,  0, 0, -1)
  A <- Matrix::Matrix(matrix(A_data, nrow = 4, ncol = 3), sparse = TRUE)
  b <- c(1, 0, 0, 0)
  c <- c(1, 2, 3)
  cone_dict <- list(z = 1L, l = 3L)
  ctrl <- list(tol_gap_abs = 1e-12, tol_gap_rel = 1e-12, tol_feas = 1e-12)

  res <- do.call(diffcp::solve_and_derivative,
                 c(list(A, b, c, cone_dict, mode = "lsqr"), ctrl))

  A_csc <- methods::as(A, "CsparseMatrix")
  nz <- Matrix::summary(A_csc)
  dA <- Matrix::sparseMatrix(i = nz$i, j = nz$j,
                             x = rnorm(nrow(nz)) * 1e-6,
                             dims = c(nrow(A), ncol(A)))
  db <- rnorm(length(b)) * 1e-6
  dc <- rnorm(length(c)) * 1e-6

  d <- res$D(dA, db, dc)
  res_pert <- do.call(diffcp::solve_only,
                      c(list(A + dA, b + db, c + dc, cone_dict), ctrl))

  ## LSQR is an iterative solve; its convergence tolerance bounds the
  ## achievable agreement.  Use the same atol=1e-8 as Python (Python's
  ## test_clarabel.py uses mode="lsqr" by default and asserts 1e-8).
  expect_lt(max(abs(d$dx - (res_pert$x - res$x))), 1e-8)
  expect_lt(max(abs(d$dy - (res_pert$y - res$y))), 1e-8)
  expect_lt(max(abs(d$ds - (res_pert$s - res$s))), 1e-8)
})
