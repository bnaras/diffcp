## CVXPY parity: each test pins to numerics computed via Python
## diffcp's cones.pi / _proj on the same input vector. See
## scripts/diffcp_expected.py (in the diffcp R package source) for
## how to regenerate.

# -- Zero cone -----------------------------------------------------

test_that("Zero cone primal projection is zero", {
  expect_equal(diffcp:::.proj_zero(c(1, 2, 3)), c(0, 0, 0))
})

test_that("Zero cone dual projection is identity", {
  expect_equal(diffcp:::.proj_zero(c(1, 2, 3), dual = TRUE), c(1, 2, 3))
})

# -- Nonneg orthant ------------------------------------------------

test_that("Nonneg orthant projection clips negatives to 0", {
  expect_equal(diffcp:::.proj_nonneg(c(-1, 0, 2, -3)), c(0, 0, 2, 0))
})

# -- SOC -----------------------------------------------------------

test_that("SOC projection: interior of K is the identity", {
  ## ||z|| < t -> in the cone -> identity
  x <- c(2, 0.5, 0.5)
  expect_equal(diffcp:::.proj_soc_one(x), x)
})

test_that("SOC projection: interior of -K is zero", {
  ## ||z|| < -t -> in the polar -> zero
  x <- c(-2, 0.5, 0.5)
  expect_equal(diffcp:::.proj_soc_one(x), c(0, 0, 0))
})

test_that("SOC projection: transition region uses the alpha formula", {
  ## ||z|| > |t| -> alpha = (t + ||z||)/2; out = (alpha, alpha * z / ||z||).
  ## With t=1, z=(3, 4): ||z||=5, alpha=3, scale=alpha/||z||=0.6.
  ## out = (3, 1.8, 2.4).
  x <- c(1, 3, 4)
  expect_equal(diffcp:::.proj_soc_one(x), c(3, 1.8, 2.4), tolerance = 1e-10)
})

# -- pi() over a Cartesian product --------------------------------

test_that("pi() composes per-cone projections in SCS order", {
  ## A small mixed problem: 2 zero, 3 nonneg, 1 SOC of size 3.
  cone_dict <- list(z = 2, l = 3, q = 3)
  cones <- parse_cone_dict(cone_dict)
  x <- c(1, 2,            # zero
         -1, 0, 5,        # nonneg
         1, 3, 4)         # SOC
  expected <- c(0, 0,     # zero -> 0
                0, 0, 5,  # nonneg -> max(., 0)
                3, 1.8, 2.4)  # SOC interior of K?  No: ||(3,4)||=5 > 1
  expect_equal(pi(x, cones), expected, tolerance = 1e-10)
})

# -- vec_symm / unvec_symm round-trip -----------------------------

test_that("vec_symm and unvec_symm round-trip", {
  X <- matrix(c(1, 2, 3,
                2, 4, 5,
                3, 5, 6), 3, 3, byrow = TRUE)
  v <- vec_symm(X)
  expect_length(v, 6L)  # 3*(3+1)/2 = 6
  X2 <- unvec_symm(v, 3L)
  expect_equal(X2, X, tolerance = 1e-10)
})
