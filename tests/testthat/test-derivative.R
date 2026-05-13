## Derivative parity vs Python `diffcp.cone_program.solve_and_derivative`
## (mode = "dense", Clarabel forward solve). Expected (dx, dy, ds) and
## adjoint (dA, db, dc) values were generated via:
##
##     uv run python inst/python_parity/derivative.py
##
## All sparse-matrix data values are pinned in CSC (column-major)
## order so the R-side `Matrix::summary()` round-trip matches without
## any reordering.

skip_if_not_installed("clarabel")

# -- LP fixture: x in {0} x R_+^3, eq + nonneg ---------------------
test_that("D and DT match Python diffcp on the LP fixture", {
  A_data <- c(1, -1, 0, 0,
              1,  0, -1, 0,
              1,  0, 0, -1)
  A <- Matrix::Matrix(matrix(A_data, nrow = 4, ncol = 3), sparse = TRUE)
  b <- c(1, 0, 0, 0)
  c <- c(1, 2, 3)
  cone_dict <- list(z = 1L, l = 3L)

  res <- diffcp::solve_and_derivative(A, b, c, cone_dict,
                                      mode = "dense",
                                      tol_gap_abs = 1e-10,
                                      tol_gap_rel = 1e-10,
                                      tol_feas    = 1e-10)

  ## CSC nz pattern of A: col1=(1,1)(2,1) col2=(1,2)(3,2) col3=(1,3)(4,3)
  dA <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 3, 1, 4),
    j = c(1, 1, 2, 2, 3, 3),
    x = c(0.000496714153011, 0.00152302985641, -0.000138264301171,
         -0.000234153374723, 0.000647688538101, -0.000234136956949),
    dims = c(4, 3)
  )
  db <- c(0.001, 0.002, 0.003, 0.004)
  dc <- c(0.001, 0.002, 0.003)

  d <- res$D(dA, db, dc)
  expect_equal(d$dx, c(0.00750328584711, -0.00300000000009, -0.00400000000001),
               tolerance = 1e-10)
  expect_equal(d$dy, c(-0.000503285846809, -4.3247428433e-14,
                       0.00140082507945, 0.00138075170098),
               tolerance = 1e-10)
  expect_equal(d$ds, c(0, 0.00798025599067,
                       -1.63432926233e-13, -9.98264453618e-14),
               tolerance = 1e-10)

  ## Adjoint derivative.
  dxa <- c(0.257399925345, -0.908481432781, -0.378503106059)
  dya <- c(-0.534915598776, 0.858073346072, -0.413009982315, 0.498188584487)
  dsa <- c(2.01019924757, 1.26286154452, -0.439214856868, -0.346437892989)
  da <- res$DT(dxa, dya, dsa)
  dA_csc <- methods::as(da$dA, "CsparseMatrix")
  ## CSC `@x` is already column-major.
  expect_equal(dA_csc@x,
               c(-1.96999846653, -1.26286154449,
                  0.413009982335, -0.413009982378,
                 -0.498188584491,  0.996377168926),
               tolerance = 1e-10)
  expect_equal(da$db,
               c(1.5202614699, 1.26286154451, 2.42874290267, 1.89876457593),
               tolerance = 1e-10)
  expect_equal(da$dc,
               c(0.44973699665, -0.413009982348, 0.49818858447),
               tolerance = 1e-10)
})

# -- SOC fixture: minimize t s.t. ||x||<=t, sum(x)=1 ---------------
test_that("D and DT match Python diffcp on the SOC fixture", {
  A_data <- c(0, -1, 0, 0, 0,
              1,  0, -1, 0, 0,
              1,  0, 0, -1, 0,
              1,  0, 0, 0, -1)
  A <- Matrix::Matrix(matrix(A_data, nrow = 5, ncol = 4), sparse = TRUE)
  b <- c(1, 0, 0, 0, 0)
  c <- c(1, 0, 0, 0)
  cone_dict <- list(z = 1L, q = 4L)

  res <- diffcp::solve_and_derivative(A, b, c, cone_dict,
                                      mode = "dense",
                                      tol_gap_abs = 1e-10,
                                      tol_gap_rel = 1e-10,
                                      tol_feas    = 1e-10)

  ## CSC nz pattern of A:
  ##   col1: (2,1); col2: (1,2)(3,2); col3: (1,3)(4,3); col4: (1,4)(5,4)
  dA <- Matrix::sparseMatrix(
    i = c(2, 1, 3, 1, 4, 1, 5),
    j = c(1, 2, 2, 3, 3, 4, 4),
    x = c(0.00152302985641, 0.000496714153011, -0.000234153374723,
         -0.000138264301171, -0.000234136956949, 0.000647688538101,
          0.00157921281551),
    dims = c(5, 4)
  )
  db <- c(0.001, 0.002, 0.003, 0.004, 0.005)
  dc <- c(0.001, 0.002, 0.003, 0.004)

  d <- res$D(dA, db, dc)
  expect_equal(d$dx,
               c(0.00597744664147, 0.00144969476511,
                -0.000339304043643, -0.000445770184804),
               tolerance = 1e-10)
  expect_equal(d$dy,
               c(-0.00404924341143, 0.00252302985636,
                -0.00220083294743, -0.000834237444808, -0.00133494550732),
               tolerance = 1e-10)
  expect_equal(d$ds,
               c(0, 0.00709812494456, 0.00452774588999,
                 0.00373874160865, 0.00402782554334),
               tolerance = 1e-10)

  ## Adjoint derivative.
  dxa <- c(0.257399925345, -0.908481432781, -0.378503106059, -0.534915598776)
  dya <- c(0.858073346072, -0.413009982315,
           0.498188584487, 2.01019924757, 1.26286154452)
  dsa <- c(-0.439214856868, -0.346437892989,
           0.455319659557, -1.66866270701, -0.86208549501)
  da <- res$DT(dxa, dya, dsa)
  dA_csc <- methods::as(da$dA, "CsparseMatrix")
  expect_equal(dA_csc@x,
               c(-2.93714078599,
                  1.33543874044,  1.03261159618,
                 -0.068855702804, -0.195023404827,
                  0.579341494297, 0.401036294707),
               tolerance = 1e-10)
  expect_equal(da$db,
               c(-1.35051565463, -0.257399925291,
                -0.442034221847, -0.972012548569, -0.815600055853),
               tolerance = 1e-10)
  expect_equal(da$dc,
               c(-3.08575070228, -1.53332717217,
                 0.89898215244, -0.223728326325),
               tolerance = 1e-10)
})

test_that("mode='lsqr' agrees with mode='dense' on D and DT (LP)", {
  A_data <- c(1, -1, 0, 0,
              1,  0, -1, 0,
              1,  0, 0, -1)
  A <- Matrix::Matrix(matrix(A_data, nrow = 4, ncol = 3), sparse = TRUE)
  b <- c(1, 0, 0, 0)
  c <- c(1, 2, 3)
  cone_dict <- list(z = 1L, l = 3L)
  ctrl <- list(tol_gap_abs = 1e-10, tol_gap_rel = 1e-10, tol_feas = 1e-10)

  res_d <- do.call(diffcp::solve_and_derivative,
                   c(list(A, b, c, cone_dict, mode = "dense"), ctrl))
  res_l <- do.call(diffcp::solve_and_derivative,
                   c(list(A, b, c, cone_dict, mode = "lsqr"),  ctrl))

  dA <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 3, 1, 4),
    j = c(1, 1, 2, 2, 3, 3),
    x = c(0.000496714153011, 0.00152302985641, -0.000138264301171,
         -0.000234153374723, 0.000647688538101, -0.000234136956949),
    dims = c(4, 3)
  )
  db <- c(0.001, 0.002, 0.003, 0.004)
  dc <- c(0.001, 0.002, 0.003)
  dD_dense <- res_d$D(dA, db, dc)
  dD_lsqr  <- res_l$D(dA, db, dc)
  expect_equal(dD_lsqr$dx, dD_dense$dx, tolerance = 1e-6)
  expect_equal(dD_lsqr$dy, dD_dense$dy, tolerance = 1e-6)
  expect_equal(dD_lsqr$ds, dD_dense$ds, tolerance = 1e-6)

  dxa <- c(0.257399925345, -0.908481432781, -0.378503106059)
  dya <- c(-0.534915598776, 0.858073346072, -0.413009982315, 0.498188584487)
  dsa <- c(2.01019924757, 1.26286154452, -0.439214856868, -0.346437892989)
  dT_dense <- res_d$DT(dxa, dya, dsa)
  dT_lsqr  <- res_l$DT(dxa, dya, dsa)
  expect_equal(as.vector(methods::as(dT_lsqr$dA, "CsparseMatrix")@x),
               as.vector(methods::as(dT_dense$dA, "CsparseMatrix")@x),
               tolerance = 1e-6)
  expect_equal(dT_lsqr$db, dT_dense$db, tolerance = 1e-6)
  expect_equal(dT_lsqr$dc, dT_dense$dc, tolerance = 1e-6)
})
