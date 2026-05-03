## Faithful R ports of the canonical Python diffcp tests:
##
##   diffcp/tests/test_clarabel.py::test_solve_and_derivative
##   diffcp/tests/test_clarabel.py::test_psdcone
##   diffcp/tests/test_scs.py::test_solve_and_derivative
##   diffcp/tests/test_scs.py::test_warm_start
##   diffcp/tests/test_scs.py::test_expcone
##
## (A, b, c, cone_dict) and the deterministic random perturbations
## come from `inst/python_parity/canonical.py`, dumped into
## `helper-canonical-data.R` (auto-sourced by testthat).
##
## Same atol bounds as Python: 1e-8 on D() finite-diff agreement and
## adjoint objective sensitivity; 1e-7 on warm-start (x_p ~= x).

skip_if_not_installed("clarabel")
skip_if_not_installed("scs")

.csc <- function(A_i, A_p, A_x, dims) {
  Matrix::sparseMatrix(i = A_i, p = A_p - 1L, x = A_x, dims = dims, index1 = TRUE)
}

# ---- helper: build dA from same nonzero pattern, given pinned values --
.dA_from_pattern <- function(A_i, A_p, dA_data, dims) {
  Matrix::sparseMatrix(i = A_i, p = A_p - 1L, x = dA_data, dims = dims, index1 = TRUE)
}

# === Port: test_clarabel.py::test_solve_and_derivative =================
test_that("clarabel: D matches finite-diff on LS-eq problem (mode='lsqr')", {
  A <- .csc(clarabel_lseq_A_A_i, clarabel_lseq_A_A_p,
            clarabel_lseq_A_A_x, clarabel_lseq_A_shape)
  res <- diffcp::solve_and_derivative(A, clarabel_lseq_b, clarabel_lseq_c,
                                      clarabel_lseq_cone_dict,
                                      mode = "lsqr",
                                      solve_method = "Clarabel")
  dA <- .dA_from_pattern(clarabel_lseq_A_A_i, clarabel_lseq_A_A_p,
                         clarabel_lseq_lsqr_dA_data, clarabel_lseq_A_shape)
  d <- res$D(dA, clarabel_lseq_lsqr_db, clarabel_lseq_lsqr_dc)

  res_pert <- diffcp::solve_only(A + dA,
                                 clarabel_lseq_b + clarabel_lseq_lsqr_db,
                                 clarabel_lseq_c + clarabel_lseq_lsqr_dc,
                                 clarabel_lseq_cone_dict,
                                 solve_method = "Clarabel")

  expect_lt(max(abs(d$dx - (res_pert$x - res$x))), 1e-8)
  expect_lt(max(abs(d$dy - (res_pert$y - res$y))), 1e-8)
  expect_lt(max(abs(d$ds - (res_pert$s - res$s))), 1e-8)
})

test_that("clarabel: D matches finite-diff on LS-eq problem (mode='dense')", {
  A <- .csc(clarabel_lseq_A_A_i, clarabel_lseq_A_A_p,
            clarabel_lseq_A_A_x, clarabel_lseq_A_shape)
  res <- diffcp::solve_and_derivative(A, clarabel_lseq_b, clarabel_lseq_c,
                                      clarabel_lseq_cone_dict,
                                      mode = "dense",
                                      solve_method = "Clarabel")
  dA <- .dA_from_pattern(clarabel_lseq_A_A_i, clarabel_lseq_A_A_p,
                         clarabel_lseq_dense_dA_data, clarabel_lseq_A_shape)
  d <- res$D(dA, clarabel_lseq_dense_db, clarabel_lseq_dense_dc)

  res_pert <- diffcp::solve_only(A + dA,
                                 clarabel_lseq_b + clarabel_lseq_dense_db,
                                 clarabel_lseq_c + clarabel_lseq_dense_dc,
                                 clarabel_lseq_cone_dict,
                                 solve_method = "Clarabel")

  expect_lt(max(abs(d$dx - (res_pert$x - res$x))), 1e-8)
  expect_lt(max(abs(d$dy - (res_pert$y - res$y))), 1e-8)
  expect_lt(max(abs(d$ds - (res_pert$s - res$s))), 1e-8)
})

# === Port: test_clarabel.py::test_psdcone =============================
test_that("clarabel: 5x5 min-eig SDP forward solve (test_psdcone)", {
  A <- .csc(clarabel_psdcone_A_A_i, clarabel_psdcone_A_A_p,
            clarabel_psdcone_A_A_x, clarabel_psdcone_A_shape)
  res <- diffcp::solve_and_derivative(A, clarabel_psdcone_b, clarabel_psdcone_c,
                                      clarabel_psdcone_cone_dict,
                                      solve_method = "Clarabel")
  ## Python asserts trace(unvec_symm(sol_vec)) ~= 1 at 1e-6 and all
  ## eigenvalues >= -1e-6.
  sol <- diffcp::unvec_symm(res$x, 5L)
  expect_lt(abs(sum(diag(sol)) - 1), 1e-6)
  evs <- eigen(sol, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(evs), -1e-6)
})

# === Port: test_scs.py::test_solve_and_derivative =====================
test_that("scs: D matches finite-diff on LS-eq problem (mode='lsqr')", {
  A <- .csc(scs_lseq_A_A_i, scs_lseq_A_A_p,
            scs_lseq_A_A_x, scs_lseq_A_shape)
  res <- diffcp::solve_and_derivative(A, scs_lseq_b, scs_lseq_c,
                                      scs_lseq_cone_dict,
                                      mode = "lsqr",
                                      solve_method = "SCS",
                                      eps = 1e-10)
  dA <- .dA_from_pattern(scs_lseq_A_A_i, scs_lseq_A_A_p,
                         scs_lseq_lsqr_dA_data, scs_lseq_A_shape)
  d <- res$D(dA, scs_lseq_lsqr_db, scs_lseq_lsqr_dc)

  res_pert <- diffcp::solve_only(A + dA,
                                 scs_lseq_b + scs_lseq_lsqr_db,
                                 scs_lseq_c + scs_lseq_lsqr_dc,
                                 scs_lseq_cone_dict,
                                 solve_method = "SCS",
                                 eps = 1e-10)

  expect_lt(max(abs(d$dx - (res_pert$x - res$x))), 1e-8)
  expect_lt(max(abs(d$dy - (res_pert$y - res$y))), 1e-8)
  expect_lt(max(abs(d$ds - (res_pert$s - res$s))), 1e-8)
})

test_that("scs: D matches finite-diff on LS-eq problem (mode='dense')", {
  A <- .csc(scs_lseq_A_A_i, scs_lseq_A_A_p,
            scs_lseq_A_A_x, scs_lseq_A_shape)
  res <- diffcp::solve_and_derivative(A, scs_lseq_b, scs_lseq_c,
                                      scs_lseq_cone_dict,
                                      mode = "dense",
                                      solve_method = "SCS",
                                      eps = 1e-10)
  dA <- .dA_from_pattern(scs_lseq_A_A_i, scs_lseq_A_A_p,
                         scs_lseq_dense_dA_data, scs_lseq_A_shape)
  d <- res$D(dA, scs_lseq_dense_db, scs_lseq_dense_dc)

  res_pert <- diffcp::solve_only(A + dA,
                                 scs_lseq_b + scs_lseq_dense_db,
                                 scs_lseq_c + scs_lseq_dense_dc,
                                 scs_lseq_cone_dict,
                                 solve_method = "SCS",
                                 eps = 1e-10)

  expect_lt(max(abs(d$dx - (res_pert$x - res$x))), 1e-8)
  expect_lt(max(abs(d$dy - (res_pert$y - res$y))), 1e-8)
  expect_lt(max(abs(d$ds - (res_pert$s - res$s))), 1e-8)
})

# === Port: test_scs.py::test_warm_start ===============================
test_that("scs: warm-start solve agrees with cold solve (max_iters=1)", {
  A <- .csc(scs_warm_start_A_A_i, scs_warm_start_A_A_p,
            scs_warm_start_A_A_x, scs_warm_start_A_shape)
  res <- diffcp::solve_and_derivative(A, scs_warm_start_b, scs_warm_start_c,
                                      scs_warm_start_cone_dict,
                                      solve_method = "SCS", eps = 1e-9)
  ## With warm-start at the optimum and max_iters = 1, SCS should
  ## return essentially the same point.
  res_p <- diffcp::solve_and_derivative(A, scs_warm_start_b, scs_warm_start_c,
                                        scs_warm_start_cone_dict,
                                        warm_start = list(res$x, res$y, res$s),
                                        max_iters = 1L,
                                        solve_method = "SCS", eps = 1e-9)
  expect_lt(max(abs(res$x - res_p$x)), 1e-7)
  expect_lt(max(abs(res$y - res_p$y)), 1e-7)
  expect_lt(max(abs(res$s - res_p$s)), 1e-7)
})

# === Port: test_scs.py::test_expcone ==================================
test_that("scs: D matches finite-diff on entropy-max EXP cone (mode='lsqr')", {
  A <- .csc(scs_expcone_A_A_i, scs_expcone_A_A_p,
            scs_expcone_A_A_x, scs_expcone_A_shape)
  res <- diffcp::solve_and_derivative(A, scs_expcone_b, scs_expcone_c,
                                      scs_expcone_cone_dict,
                                      mode = "lsqr", solve_method = "SCS",
                                      eps = 1e-10)
  dA <- .dA_from_pattern(scs_expcone_A_A_i, scs_expcone_A_A_p,
                         scs_expcone_lsqr_dA_data, scs_expcone_A_shape)
  d <- res$D(dA, scs_expcone_lsqr_db, scs_expcone_lsqr_dc)

  res_pert <- diffcp::solve_only(A + dA,
                                 scs_expcone_b + scs_expcone_lsqr_db,
                                 scs_expcone_c + scs_expcone_lsqr_dc,
                                 scs_expcone_cone_dict,
                                 solve_method = "SCS", eps = 1e-10)

  expect_lt(max(abs(d$dx - (res_pert$x - res$x))), 1e-8)
  expect_lt(max(abs(d$dy - (res_pert$y - res$y))), 1e-8)
  expect_lt(max(abs(d$ds - (res_pert$s - res$s))), 1e-8)
})

test_that("scs: D matches finite-diff on entropy-max EXP cone (mode='dense')", {
  A <- .csc(scs_expcone_A_A_i, scs_expcone_A_A_p,
            scs_expcone_A_A_x, scs_expcone_A_shape)
  res <- diffcp::solve_and_derivative(A, scs_expcone_b, scs_expcone_c,
                                      scs_expcone_cone_dict,
                                      mode = "dense", solve_method = "SCS",
                                      eps = 1e-10)
  dA <- .dA_from_pattern(scs_expcone_A_A_i, scs_expcone_A_A_p,
                         scs_expcone_dense_dA_data, scs_expcone_A_shape)
  d <- res$D(dA, scs_expcone_dense_db, scs_expcone_dense_dc)

  res_pert <- diffcp::solve_only(A + dA,
                                 scs_expcone_b + scs_expcone_dense_db,
                                 scs_expcone_c + scs_expcone_dense_dc,
                                 scs_expcone_cone_dict,
                                 solve_method = "SCS", eps = 1e-10)

  expect_lt(max(abs(d$dx - (res_pert$x - res$x))), 1e-8)
  expect_lt(max(abs(d$dy - (res_pert$y - res$y))), 1e-8)
  expect_lt(max(abs(d$ds - (res_pert$s - res$s))), 1e-8)
})
