## Faithful R ports of diffcp/tests/test_cone_prog_diff.py.
##
## Tests grouped by topic (matching the order in the Python file):
##   * cone-dimension utilities (vec_psd_dim, psd_dim)
##   * EXP-cone membership (in_exp, in_exp_dual)
##   * vec_symm / unvec_symm round-trips
##   * per-cone projection sanity (_proj_zero/pos/soc/psd/exp)
##   * per-cone Jacobian sanity (dproj_*) via FD
##   * pi() and dpi() composition over multi-cone product
##   * lsqr vs lstsq sanity
##   * infeasible-problem error
##   * dprojection_exp analytical sensitivity
##
## Same atol bounds as Python.  Fixtures requiring cvxpy reference
## values come from inst/python_parity/cone_diff.py via
## helper-cone-diff-data.R (auto-sourced by testthat).

skip_if_not_installed("clarabel")
skip_if_not_installed("scs")

# ---- helper -----------------------------------------------------
.csc <- function(A_i, A_p, A_x, dims) {
  Matrix::sparseMatrix(i = A_i, p = A_p - 1L, x = A_x, dims = dims, index1 = TRUE)
}

# === test_vec_psd_dim ============================================
test_that("vec_psd_dim(10) == 55", {
  ## Mirrors Python `cone_lib.vec_psd_dim(10) == 10*11/2`.
  expect_equal(diffcp:::.vec_psd_dim(10L), 55L)
})

# === test_psd_dim ================================================
test_that("psd_dim(vec_psd_dim(4096)) == 4096", {
  expect_equal(diffcp:::.psd_dim(diffcp:::.vec_psd_dim(4096L)), 4096L)
})

# === test_in_exp =================================================
test_that("in_exp matches Python", {
  expect_true(diffcp:::cpp_in_exp(c(0, 0, 1)))
  expect_true(diffcp:::cpp_in_exp(c(-1, 0, 0)))
  expect_true(diffcp:::cpp_in_exp(c(1, 1, 5)))
  expect_false(diffcp:::cpp_in_exp(c(1, 0, 0)))
  expect_false(diffcp:::cpp_in_exp(c(-1, -1, 1)))
  expect_false(diffcp:::cpp_in_exp(c(-1, 0, -1)))
})

# === test_in_exp_dual ============================================
test_that("in_exp_dual matches Python", {
  expect_true(diffcp:::cpp_in_exp_dual(c(0, 1, 1)))
  expect_true(diffcp:::cpp_in_exp_dual(c(-1, 1, 5)))
  expect_false(diffcp:::cpp_in_exp_dual(c(0, -1, 1)))
  expect_false(diffcp:::cpp_in_exp_dual(c(0, 1, -1)))
})

# === test_unvec_symm =============================================
test_that("unvec_symm(vec_symm(X)) == X for symmetric X", {
  set.seed(0)
  n <- 5L
  x <- matrix(rnorm(n * n), n, n)
  x <- x + t(x)
  expect_equal(diffcp::unvec_symm(diffcp::vec_symm(x), n), x,
               tolerance = 1e-12)
})

# === test_vec_symm ===============================================
test_that("vec_symm(unvec_symm(x)) == x", {
  set.seed(0)
  n <- 5L
  m <- (n * (n + 1L)) %/% 2L
  x <- rnorm(m)
  expect_equal(diffcp::vec_symm(diffcp::unvec_symm(x, n)), x,
               tolerance = 1e-12)
})

# === test_proj_zero ==============================================
test_that("_proj_zero: dual=TRUE returns x; dual=FALSE returns 0", {
  set.seed(0)
  n <- 100L
  for (k in 1:10) {
    x <- rnorm(n)
    cones <- list(list(name = "z", sizes = n))
    expect_equal(diffcp::pi(x, cones, dual = TRUE),  x,         tolerance = 1e-15)
    expect_equal(diffcp::pi(x, cones, dual = FALSE), numeric(n),tolerance = 1e-15)
  }
})

# === test_proj_pos ===============================================
test_that("_proj_pos: returns max(x, 0); self-dual", {
  set.seed(0)
  n <- 100L
  cones <- list(list(name = "l", sizes = n))
  for (k in 1:15) {
    x <- rnorm(n)
    p <- diffcp::pi(x, cones, dual = FALSE)
    expect_equal(p, pmax(x, 0), tolerance = 1e-15)
    expect_equal(p, diffcp::pi(x, cones, dual = TRUE), tolerance = 1e-15)
  }
})

# === test_proj_soc ===============================================
test_that("_proj_soc matches Python diffcp reference (15 random vectors)", {
  ## Pinned reference from inst/python_parity/cone_diff.py:emit_proj_soc.
  n <- 100L
  cones <- list(list(name = "q", sizes = n))
  for (k in 1:15) {
    x <- proj_soc_x[((k - 1L) * n + 1L):(k * n)]
    p_ref <- proj_soc_p[((k - 1L) * n + 1L):(k * n)]
    p <- diffcp::pi(x, cones, dual = FALSE)
    expect_lt(max(abs(p - p_ref)), 1e-12)
    ## SOC is self-dual.
    expect_lt(max(abs(p - diffcp::pi(x, cones, dual = TRUE))), 1e-15)
  }
})

# === test_proj_psd ===============================================
test_that("_proj_psd matches Python diffcp reference (15 random matrices)", {
  n <- 10L
  m <- (n * (n + 1L)) %/% 2L
  cones <- list(list(name = "s", sizes = n))
  for (k in 1:15) {
    x <- proj_psd_x[((k - 1L) * m + 1L):(k * m)]
    p_ref <- proj_psd_p[((k - 1L) * m + 1L):(k * m)]
    p <- diffcp::pi(x, cones, dual = FALSE)
    expect_lt(max(abs(p - p_ref)), 1e-10)
    ## PSD is self-dual.
    expect_lt(max(abs(p - diffcp::pi(x, cones, dual = TRUE))), 1e-12)
  }
})

# === test_proj_exp ===============================================
test_that("_proj_exp matches Python diffcp reference (15 random vectors)", {
  ## 9-dim = 3 EXP cones.
  cones_exp <- list(list(name = "ep", sizes = 3L))
  cones_exp_dual <- list(list(name = "ed", sizes = 3L))
  for (k in 1:15) {
    x <- proj_exp_x[((k - 1L) * 9L + 1L):(k * 9L)]
    p_ref <- proj_exp_p[((k - 1L) * 9L + 1L):(k * 9L)]
    p_dual_ref <- proj_exp_p_dual[((k - 1L) * 9L + 1L):(k * 9L)]
    expect_lt(max(abs(diffcp::pi(x, cones_exp, dual = FALSE) - p_ref)), 1e-10)
    expect_lt(max(abs(diffcp::pi(x, cones_exp_dual, dual = FALSE) - p_dual_ref)),
              1e-10)
  }
})

# === Helper for dproj_* finite-difference test (mirrors _test_dproj) -
.test_dproj <- function(cone, dual, n, tol = 1e-8) {
  x <- rnorm(n)
  dx <- 1e-6 * rnorm(n)
  cones <- list(cone)
  proj_x <- diffcp::pi(x,      cones, dual = dual)
  proj_z <- diffcp::pi(x + dx, cones, dual = dual)
  Dpi_dense <- diffcp:::cpp_dprojection_dense(x, cones, dual)
  expect_lt(max(abs(Dpi_dense %*% dx - (proj_z - proj_x))), tol)
}

# === test_dproj_zero =============================================
test_that("dproj_zero matches FD (10 random)", {
  set.seed(0)
  cone <- list(name = "z", sizes = 55L)
  for (k in 1:10) .test_dproj(cone, dual = TRUE,  n = 55L)
  for (k in 1:10) .test_dproj(cone, dual = FALSE, n = 55L)
})

# === test_dproj_pos ==============================================
test_that("dproj_pos matches FD (10 random)", {
  set.seed(0)
  cone <- list(name = "l", sizes = 55L)
  for (k in 1:10) .test_dproj(cone, dual = TRUE,  n = 55L)
  for (k in 1:10) .test_dproj(cone, dual = FALSE, n = 55L)
})

# === test_dproj_soc ==============================================
test_that("dproj_soc matches FD (10 random)", {
  set.seed(0)
  cone <- list(name = "q", sizes = 55L)
  for (k in 1:10) .test_dproj(cone, dual = TRUE,  n = 55L)
  for (k in 1:10) .test_dproj(cone, dual = FALSE, n = 55L)
})

# === test_dproj_psd ==============================================
test_that("dproj_psd matches FD (10 random)", {
  set.seed(0)
  ## n=55 packed -> matrix dim 10
  cone <- list(name = "s", sizes = diffcp:::.psd_dim(55L))
  for (k in 1:10) .test_dproj(cone, dual = TRUE,  n = 55L)
  for (k in 1:10) .test_dproj(cone, dual = FALSE, n = 55L)
})

# === test_dproj_exp ==============================================
test_that("dproj_exp matches FD (10 random, tol=1e-5)", {
  set.seed(0)
  cone <- list(name = "ep", sizes = 18L)
  for (k in 1:10) .test_dproj(cone, dual = TRUE,  n = 54L, tol = 1e-5)
  for (k in 1:10) .test_dproj(cone, dual = FALSE, n = 54L, tol = 1e-5)
})

# === test_dproj_exp_dual =========================================
test_that("dproj_exp_dual matches FD (10 random, tol=1e-5)", {
  set.seed(0)
  cone <- list(name = "ed", sizes = 18L)
  for (k in 1:10) .test_dproj(cone, dual = TRUE,  n = 54L, tol = 1e-5)
  for (k in 1:10) .test_dproj(cone, dual = FALSE, n = 54L, tol = 1e-5)
})

# === test_pi (multi-cone composition) ============================
test_that("pi composes per-block _proj across multi-cone product", {
  set.seed(0)
  for (rep in 1:10) {
    zero_dim <- sample.int(9L, 1L)
    pos_dim  <- sample.int(9L, 1L)
    soc_dim  <- sample.int(9L, sample.int(9L, 1L), replace = TRUE)
    psd_dim  <- sample.int(9L, sample.int(9L, 1L), replace = TRUE)
    exp_dim  <- sample.int(15L, 1L) + 2L  # 3..17

    cones <- list(
      list(name = "z",  sizes = zero_dim),
      list(name = "l",  sizes = pos_dim),
      list(name = "q",  sizes = as.integer(soc_dim)),
      list(name = "s",  sizes = as.integer(psd_dim)),
      list(name = "ep", sizes = exp_dim),
      list(name = "ed", sizes = exp_dim)
    )
    size <- zero_dim + pos_dim + sum(soc_dim) +
      sum(sapply(psd_dim, function(d) (d * (d + 1L)) %/% 2L)) +
      2L * 3L * exp_dim
    x <- rnorm(size)

    for (dual in c(FALSE, TRUE)) {
      proj <- diffcp::pi(x, cones, dual = dual)

      offset <- 0L
      ## zero block
      block_z <- diffcp::pi(x[(offset + 1L):(offset + zero_dim)],
                            list(list(name = "z", sizes = zero_dim)),
                            dual = dual)
      expect_lt(max(abs(proj[(offset + 1L):(offset + zero_dim)] - block_z)), 1e-12)
      offset <- offset + zero_dim

      ## pos block
      block_l <- diffcp::pi(x[(offset + 1L):(offset + pos_dim)],
                            list(list(name = "l", sizes = pos_dim)),
                            dual = dual)
      expect_lt(max(abs(proj[(offset + 1L):(offset + pos_dim)] - block_l)), 1e-12)
      offset <- offset + pos_dim

      ## soc blocks
      for (d in soc_dim) {
        block_q <- diffcp::pi(x[(offset + 1L):(offset + d)],
                              list(list(name = "q", sizes = as.integer(d))),
                              dual = dual)
        expect_lt(max(abs(proj[(offset + 1L):(offset + d)] - block_q)), 1e-12)
        offset <- offset + d
      }

      ## psd blocks
      for (d in psd_dim) {
        m <- (d * (d + 1L)) %/% 2L
        block_s <- diffcp::pi(x[(offset + 1L):(offset + m)],
                              list(list(name = "s", sizes = as.integer(d))),
                              dual = dual)
        expect_lt(max(abs(proj[(offset + 1L):(offset + m)] - block_s)), 1e-10)
        offset <- offset + m
      }

      ## exp block
      block_ep <- diffcp::pi(x[(offset + 1L):(offset + 3L * exp_dim)],
                             list(list(name = "ep", sizes = exp_dim)),
                             dual = dual)
      expect_lt(max(abs(proj[(offset + 1L):(offset + 3L * exp_dim)] - block_ep)),
                1e-10)
      offset <- offset + 3L * exp_dim

      ## exp_dual block (everything that's left)
      block_ed <- diffcp::pi(x[(offset + 1L):length(x)],
                             list(list(name = "ed", sizes = exp_dim)),
                             dual = dual)
      expect_lt(max(abs(proj[(offset + 1L):length(x)] - block_ed)), 1e-10)
    }
  }
})

# === test_dpi (multi-cone Jacobian composition) ==================
test_that("D pi(x) ~= (pi(x+dx) - pi(x)) / FD on multi-cone (10 random)", {
  set.seed(0)
  for (rep in 1:10) {
    zero_dim <- sample.int(9L, 1L)
    pos_dim  <- sample.int(9L, 1L)
    soc_dim  <- sample.int(9L, sample.int(9L, 1L), replace = TRUE)
    psd_dim  <- sample.int(9L, sample.int(9L, 1L), replace = TRUE)
    exp_dim  <- sample.int(15L, 1L) + 2L

    cones <- list(
      list(name = "z",  sizes = zero_dim),
      list(name = "l",  sizes = pos_dim),
      list(name = "q",  sizes = as.integer(soc_dim)),
      list(name = "s",  sizes = as.integer(psd_dim)),
      list(name = "ep", sizes = exp_dim),
      list(name = "ed", sizes = exp_dim)
    )
    size <- zero_dim + pos_dim + sum(soc_dim) +
      sum(sapply(psd_dim, function(d) (d * (d + 1L)) %/% 2L)) +
      2L * 3L * exp_dim
    x <- rnorm(size)

    for (dual in c(FALSE, TRUE)) {
      proj_x <- diffcp::pi(x, cones, dual = dual)
      dx <- 1e-7 * rnorm(size)
      z <- diffcp::pi(x + dx, cones, dual = dual)
      Dpi_dense <- diffcp:::cpp_dprojection_dense(x, cones, dual)
      expect_lt(max(abs(Dpi_dense %*% dx - (z - proj_x))), 1e-6)
    }
  }
})

# === test_lsqr ===================================================
test_that("LSQR on random A, b matches base lstsq (qr.solve)", {
  set.seed(0)
  A <- matrix(rnorm(20 * 10), 20, 10)
  b <- rnorm(20)
  A_csc <- methods::as(Matrix::Matrix(A, sparse = TRUE), "CsparseMatrix")
  out <- diffcp:::cpp_lsqr_sparse(A_csc, b)
  svx <- qr.solve(A, b)
  expect_lt(max(abs(svx - out$solution)), 1e-8)
})

# === test_infeasible =============================================
test_that("solve_only errors on a primal-infeasible 1D LP (SCS)", {
  ## c = [1], A = [[1], [1]], b = [1, -1] -> x must be both 1 and -1.
  ## Python uses cone_dims = {"f": 2} (free / equality cone).
  A <- methods::as(Matrix::Matrix(matrix(c(1, 1), nrow = 2L), sparse = TRUE),
                   "CsparseMatrix")
  b <- c(1.0, -1.0)
  c <- c(1.0)
  expect_error(
    diffcp::solve_only(A, b, c, list(z = 2L), solve_method = "SCS"),
    "SCS returned"
  )
})

# === test_dprojection_exp (analytical sensitivity) ===============
test_that("D matches analytical sensitivity for entropy-parametric problem", {
  ## CVXPY problem (one variable + parameter lam):
  ##    maximize x + lam*(log(1+x) + log(1-x))
  ## At lam = 1, sensitivity dx/dlam = -1 + lam/sqrt(lam^2 + 1) = -1 + 1/sqrt(2).
  A <- .csc(dprojexp_A_A_i, dprojexp_A_A_p, dprojexp_A_A_x, dprojexp_A_shape)

  res <- diffcp::solve_and_derivative(
    A, dprojexp_b, dprojexp_c, dprojexp_cone_dict,
    solve_method = "Clarabel", mode = "dense",
    tol_gap_abs = 1e-13, tol_gap_rel = 1e-13,
    tol_feas = 1e-13, tol_ktratio = 1e-13)

  dlam <- 1e-6
  dA <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                             dims = dim(A))
  db <- numeric(length(dprojexp_b))
  dc <- numeric(length(dprojexp_c))
  ## Python comment (cone_prog_diff.py:388-389): the minus sign stems
  ## from the fact that c = [-1, -1, -1].
  dc[2] <- -dlam
  dc[3] <- -dlam
  d <- res$D(dA, db, dc)

  expect_lt(abs(dprojexp_analytical - d$dx[1] / dlam), 1e-6)
})
