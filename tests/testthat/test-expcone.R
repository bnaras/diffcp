## EXP-cone test, ported from
##   diffcp/tests/test_clarabel.py::test_expcone
##
## Builds the entropy-maximization problem
##   min  -sum(entr(y))
##   s.t.  sum(y) = 1,  y > 0   (n = 10)
## canonicalizes to SCS form (one z cone + 10 EXP cones), solves via
## Clarabel through diffcp in modes "dense" and "lsqr" (Python also
## checks "lsmr" which we don't expose), and asserts:
##   * D(dA, db, dc) ~= (x_pert - x, y_pert - y, s_pert - s)  at atol=1e-8
## using the exact perturbations Python generates.  Same atol bound as
## test_clarabel.py:121-123.
##
## (A, b, c, cone_dict) and the perturbations come from
## inst/python_parity/expcone.py.

skip_if_not_installed("clarabel")

.csc <- function(A_i, A_p, A_x, dims) {
  Matrix::sparseMatrix(i = A_i, p = A_p - 1L, x = A_x, dims = dims, index1 = TRUE)
}

# Pinned canonical (A, b, c) for the entropy problem (n = 10) -----
A_i <- c(2, 5, 8, 11, 14, 17, 20, 23, 26, 29,
         1, 3, 1, 6, 1, 9, 1, 12, 1, 15,
         1, 18, 1, 21, 1, 24, 1, 27, 1, 30)
A_p <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
         13, 15, 17, 19, 21, 23, 25, 27, 29, 31)
A_x <- c(rep(-1, 10),
         1, -1, 1, -1, 1, -1, 1, -1, 1, -1,
         1, -1, 1, -1, 1, -1, 1, -1, 1, -1)
A <- .csc(A_i, A_p, A_x, c(31, 20))
b <- c(1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
       1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1)
c <- c(rep(-1, 10), rep(0, 10))
cone_dict <- list(z = 1L, ep = 10L)

ctrl <- list(tol_gap_abs = 1e-13, tol_gap_rel = 1e-13,
             tol_feas = 1e-13, tol_ktratio = 1e-13)

# Pinned perturbations (CSC-order dA values + db + dc) ------------
dA_data <- c(1.44043571160878e-07, 7.610377251469934e-07, 4.438632327454257e-07,
             1.494079073157606e-06, 3.130677016509014e-07, -2.552989815834078e-06,
             8.644361988595055e-07, 2.269754623987607e-06, 4.575851730144606e-08,
             1.532779214358458e-06, 1.764052345967664e-06, 1.454273506962975e-06,
             4.001572083672233e-07, 1.216750164928284e-07, 9.787379841057391e-07,
             3.336743273742668e-07, 2.240893199201458e-06, -2.051582637658009e-07,
             1.867557990149967e-06, -8.540957393017248e-07, -9.772778798764111e-07,
             6.536185954403606e-07, 9.500884175255894e-07, -7.421650204064417e-07,
             -1.513572082976979e-07, -1.454365674598765e-06,
             -1.032188517935578e-07, -1.871838500258336e-07,
             4.105985019383723e-07, 1.469358769900285e-06)
db <- c(1.549474256969163e-07, 3.781625196021735e-07, -8.877857476301127e-07,
        -1.980796468223927e-06, -3.479121493261526e-07, 1.563489691039801e-07,
        1.230290680727721e-06, 1.202379848784411e-06, -3.873268174079523e-07,
        -3.023027505753356e-07, -1.048552965067093e-06, -1.420017937178975e-06,
        -1.706270190625013e-06, 1.95077539523179e-06, -5.096521817516535e-07,
        -4.380743016111864e-07, -1.252795360049926e-06, 7.774903558319103e-07,
        -1.613897847557951e-06, -2.127402802139687e-07, -8.954665611936755e-07,
        3.86902497859262e-07, -5.10805137568873e-07, -1.180632184122412e-06,
        -2.818222833865487e-08, 4.283318705304176e-07, 6.651722238316788e-08,
        3.024718977397814e-07, -6.343220936809636e-07, -3.627411659871381e-07,
        -6.72460447775951e-07)
dc <- c(-3.595531615405413e-07, -8.13146282044454e-07, -1.726282602331677e-06,
        1.774261422537528e-07, -4.017809362082619e-07, -1.630198346966044e-06,
        4.627822555257741e-07, -9.072983643832421e-07, 5.194539579613895e-08,
        7.290905621775368e-07, 1.289829107574107e-07, 1.139400684543301e-06,
        -1.234825820353653e-06, 4.02341641177549e-07, -6.848100909403131e-07,
        -8.707971491818817e-07, -5.788496647644154e-07, -3.115525321273727e-07,
        5.616534222974543e-08, -1.165149840783356e-06)

dA <- .csc(A_i, A_p, dA_data, c(31, 20))

run_mode <- function(mode) {
  res <- do.call(diffcp::solve_and_derivative,
                 c(list(A, b, c, cone_dict,
                        solve_method = "Clarabel", mode = mode),
                   ctrl))
  d <- res$D(dA, db, dc)
  res_pert <- do.call(diffcp::solve_only,
                      c(list(A + dA, b + db, c + dc, cone_dict,
                             solve_method = "Clarabel"),
                        ctrl))
  list(res = res, d = d, res_pert = res_pert)
}

test_that("EXP cone: D matches finite-difference (mode='dense')", {
  out <- run_mode("dense")
  ## test_clarabel.py:121-123 uses atol=1e-8 for x_pert - x ~= dx etc.
  expect_lt(max(abs(out$d$dx - (out$res_pert$x - out$res$x))), 1e-8)
  expect_lt(max(abs(out$d$dy - (out$res_pert$y - out$res$y))), 1e-8)
  expect_lt(max(abs(out$d$ds - (out$res_pert$s - out$res$s))), 1e-8)
})

test_that("EXP cone: D matches finite-difference (mode='lsqr')", {
  out <- run_mode("lsqr")
  expect_lt(max(abs(out$d$dx - (out$res_pert$x - out$res$x))), 1e-8)
  expect_lt(max(abs(out$d$dy - (out$res_pert$y - out$res$y))), 1e-8)
  expect_lt(max(abs(out$d$ds - (out$res_pert$s - out$res$s))), 1e-8)
})
