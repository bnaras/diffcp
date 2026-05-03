# diffcp

[![R-CMD-check](https://github.com/bnaras/diffcp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bnaras/diffcp/actions/workflows/R-CMD-check.yaml)

`diffcp` is an R port of the Python
[`diffcp`](https://github.com/cvxgrp/diffcp) package. It computes the
derivative of the optimal solution of a convex cone program with
respect to its problem data, treating the program as an implicit
function. The derivative is exposed both as a forward operator
`D(dA, db, dc)` and an adjoint `DT(dx, dy, ds)`, and can be evaluated
in either a dense (Eigen LDLT) or matrix-free (LSQR) mode.

The implementation is a faithful port of the C++ core (cones,
LinearOperator, dprojection, M operator, LSQR) from
[cvxgrp/diffcp](https://github.com/cvxgrp/diffcp), called from R via
RcppEigen. The R-level API mirrors the Python source. Test fixtures
are pinned to Python diffcp's outputs so the two packages agree
bit-for-bit on every supported cone type.

## What it solves

Given problem data `(A, b, c)` (and optionally a quadratic objective
`P`), `diffcp` solves the primal–dual pair

```
minimise    c'x          minimise    b'y
subject to  Ax + s = b   subject to  A'y + c = 0
            s in K                   y in K*
```

where `K` is a Cartesian product of the standard cones (zero,
non-negative orthant, second-order, positive semidefinite,
exponential, exponential dual). It returns the optimal `(x, y, s)`
together with two callables:

* **`D(dA, db, dc)`** — applies the derivative of `(x, y, s)` with
  respect to `(A, b, c)` at the perturbations `(dA, db, dc)`.
* **`DT(dx, dy, ds)`** — applies the adjoint, returning the
  perturbations `(dA, db, dc)` corresponding to a desired change
  `(dx, dy, ds)` in the solution.

## Installation

```r
remotes::install_github("bnaras/diffcp")
```

`diffcp` requires a C++17 compiler and the R packages
[`Rcpp`](https://CRAN.R-project.org/package=Rcpp),
[`RcppEigen`](https://CRAN.R-project.org/package=RcppEigen),
[`Matrix`](https://CRAN.R-project.org/package=Matrix), and
[`clarabel`](https://CRAN.R-project.org/package=clarabel). The
[`scs`](https://CRAN.R-project.org/package=scs) package is an
optional alternative forward solver.

## Quick example

```r
library(diffcp)

## Tiny LP: minimise c'x s.t. 1'x = 1, x >= 0  with n = 3
A <- Matrix::sparseMatrix(
  i = c(1, 2, 1, 3, 1, 4),
  j = c(1, 1, 2, 2, 3, 3),
  x = c(1, -1, 1, -1, 1, -1),
  dims = c(4, 3))
b <- c(1, 0, 0, 0)
c <- c(1, 2, 3)
cone_dict <- list(z = 1L, l = 3L)

res <- solve_and_derivative(A, b, c, cone_dict, mode = "lsqr")
res$x   # primal optimum
res$y   # dual variable
res$s   # slack

## Apply the derivative at a perturbation in c.
dA <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                            x = numeric(0), dims = c(4, 3))
db <- numeric(4)
dc <- c(0.001, 0, 0)
res$D(dA, db, dc)   # list(dx, dy, ds)

## Apply the adjoint to dx = c (steepest descent of c'x in (A, b, c)).
res$DT(c, numeric(4), numeric(4))   # list(dA, db, dc)
```

For larger examples — including PSD and exponential cones — see the
package vignette (`vignette("diffcp", package = "diffcp")`).

## Cone support

| Cone | Tag | Forward solve | Derivative |
|------|-----|---------------|------------|
| Zero (equality) | `z` | ✓ | ✓ |
| Non-negative orthant | `l` | ✓ | ✓ |
| Second-order cone | `q` | ✓ | ✓ |
| Positive semidefinite | `s` | ✓ (Clarabel + SCS) | ✓ |
| Exponential cone | `ep` | ✓ | ✓ |
| Exponential dual cone | `ed` | ✓ | ✓ |

PSD vectorisation follows the SCS convention (lower-triangular
column-major, with off-diagonals scaled by `sqrt(2)`). When solving
PSD problems through Clarabel, `diffcp` automatically permutes the
rows of `A` and the entries of `b` from SCS order to Clarabel's
upper-triangular order, then permutes the dual `y` and slack `s` back
on return.

Quadratic objectives (`P`) are supported in `solve_only(P = P)` for
the forward solve via Clarabel or SCS native QP. The dense and LSQR
modes of `solve_and_derivative` do not support quadratic objectives;
the upstream Python package handles QP only via its `lpgd` modes,
which this R port does not yet provide.

## Differences from Python diffcp

* **No `lpgd` derivative modes**, no batch APIs, no `lsmr` mode, and
  no ECOS solver branch. The supported modes are `dense` and `lsqr`.
* **Cone projections in R, Jacobians in C++.** `pi()` and the per-cone
  projections live in R (matching Python's split between `cones.py`
  for projections and `cones.cpp` for Jacobians). All Jacobian, M,
  and LSQR machinery is in C++ via RcppEigen.
* **Eigen LDLT in dense mode.** Identical to Python's
  `(M^T M).ldlt().solve(M^T rhs)` (`cpp/src/deriv.cpp:53-57`).

## Citation

If you use `diffcp` in academic work, please cite the original paper:

```bibtex
@article{diffcp2019,
    author       = {Agrawal, A. and Barratt, S. and Boyd, S. and Busseti, E. and Moursi, W.},
    title        = {Differentiating through a cone program},
    journal      = {Journal of Applied and Numerical Optimization},
    year         = {2019},
    volume       = {1},
    number       = {2},
    pages        = {107--115},
}
```

## License

Apache License 2.0, matching the upstream Python `diffcp` license.
