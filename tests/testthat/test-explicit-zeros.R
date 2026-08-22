## Explicit zeros in `A` must survive into the pattern of `dA`.
##
## An entry of `A` that is STRUCTURALLY present but numerically zero is not the
## same thing as an absent entry: the first has a derivative, the second does
## not. Callers that differentiate with respect to a parameter driving that
## entry store it explicitly for exactly this reason -- CVXPY calls the flag
## `keep_zeros`, and it is the only reason its DIFFCP interface overrides
## `apply` (cvxpy/reductions/solvers/conic_solvers/diffcp_conif.py:72).
##
## Python diffcp preserves the pattern by marking explicit zeros as NaN,
## reading `A.nonzero()`, restoring them, and only then calling
## `eliminate_zeros()` (cone_program.py:643-653); `dA` is built on the captured
## `(rows, cols)` at :771 and :891. R needs no NaN dance -- `Matrix::summary()`
## already lists stored entries, explicit zeros included -- but the ORDER still
## matters, and this file exists because it was previously wrong: `drop0()` ran
## first and the pattern was read afterwards, so the entry vanished from `dA`
## and its gradient came back as exactly 0, with no error and no warning.

## minimize  c'x  s.t.  Ax + s = b, s >= 0, with an EXPLICIT ZERO at A[1,1].
## Feasible region x1 <= 3, x2 <= 5, and 0*x1 + x2 <= 1 -> x = (3, 1).
.ez_problem <- function(a11 = 0) {
  A <- Matrix::sparseMatrix(
    i = c(1L, 2L, 4L, 1L, 3L, 5L),
    j = c(1L, 1L, 1L, 2L, 2L, 2L),
    x = c(a11, -1,  1,  1, -1,  1),
    dims = c(5L, 2L))
  list(A = A, b = c(1, 0, 0, 3, 5), c = c(-1, -2), cone = list(z = 0L, l = 5L))
}

test_that("Matrix::summary() lists explicit zeros (the assumption behind the fix)", {
  m <- Matrix::sparseMatrix(i = c(1L, 2L), j = c(1L, 1L), x = c(0, -1),
                            dims = c(2L, 2L))
  expect_equal(nrow(Matrix::summary(m)), 2L)
  expect_equal(nrow(Matrix::summary(Matrix::drop0(m))), 1L)
})

test_that("an explicit zero in A keeps its entry in dA", {
  p <- .ez_problem(a11 = 0)
  expect_equal(length(p$A@x), 6L)         # the zero really is stored

  res <- solve_and_derivative(p$A, p$b, p$c, p$cone,
                              solve_method = "Clarabel", mode = "lsqr")
  d <- res$DT(rep(1, length(res$x)), numeric(length(res$y)), numeric(length(res$s)))
  nz <- Matrix::summary(d$dA)

  ## The pattern of dA is the pattern of A, explicit zeros included.
  expect_equal(nrow(nz), 6L)
  expect_true(any(nz$i == 1L & nz$j == 1L))

  ## And the derivative there is the real one, not a placeholder. Its value is
  ## d(sum(x))/d(A[1,1]); Python diffcp returns -3.00000001 on this problem.
  expect_equal(nz$x[nz$i == 1L & nz$j == 1L], -3, tolerance = 1e-6)
})

test_that("the derivative at an explicit zero agrees with a finite difference", {
  h <- 1e-6
  obj <- function(a11) {
    p <- .ez_problem(a11 = a11)
    r <- solve_only(Matrix::drop0(p$A), p$b, p$c, p$cone, solve_method = "Clarabel")
    sum(r$x)
  }
  fd <- (obj(h) - obj(-h)) / (2 * h)

  p <- .ez_problem(a11 = 0)
  res <- solve_and_derivative(p$A, p$b, p$c, p$cone,
                              solve_method = "Clarabel", mode = "lsqr")
  d <- res$DT(rep(1, length(res$x)), numeric(length(res$y)), numeric(length(res$s)))
  nz <- Matrix::summary(d$dA)
  analytic <- nz$x[nz$i == 1L & nz$j == 1L]

  expect_equal(analytic, fd, tolerance = 1e-4)
  ## The point of the whole exercise: it is NOT zero.
  expect_gt(abs(analytic), 1)
})

test_that("a matrix with no explicit zeros is completely unaffected", {
  ## drop0() is the identity here and summary() gives the same pattern, so the
  ## captured-pattern path must reproduce the old numbers exactly. This is what
  ## makes the change safe for every existing caller.
  p <- .ez_problem(a11 = 0)
  p$A <- Matrix::drop0(p$A)
  expect_equal(length(p$A@x), 5L)

  res <- solve_and_derivative(p$A, p$b, p$c, p$cone,
                              solve_method = "Clarabel", mode = "lsqr")
  d <- res$DT(rep(1, length(res$x)), numeric(length(res$y)), numeric(length(res$s)))
  expect_equal(nrow(Matrix::summary(d$dA)), 5L)
  expect_false(any(Matrix::summary(d$dA)$i == 1L & Matrix::summary(d$dA)$j == 1L))
})

test_that("dense mode preserves the pattern too", {
  p <- .ez_problem(a11 = 0)
  res <- solve_and_derivative(p$A, p$b, p$c, p$cone,
                              solve_method = "Clarabel", mode = "dense")
  d <- res$DT(rep(1, length(res$x)), numeric(length(res$y)), numeric(length(res$s)))
  nz <- Matrix::summary(d$dA)
  expect_equal(nrow(nz), 6L)
  expect_equal(nz$x[nz$i == 1L & nz$j == 1L], -3, tolerance = 1e-6)
})

test_that("solve_only still eliminates explicit zeros (no pattern to preserve)", {
  ## solve_only corresponds to Python cone_program.py:386, where no pattern is
  ## captured; its drop0 is correct and must stay.
  p <- .ez_problem(a11 = 0)
  r <- solve_only(p$A, p$b, p$c, p$cone, solve_method = "Clarabel")
  expect_equal(r$x, c(3, 1), tolerance = 1e-6)
})
