# Changelog

## diffcp 0.1.2

- [`solve_and_derivative()`](https://bnaras.github.io/diffcp/reference/solve_and_derivative.md)
  now reports `dA` on the sparsity pattern of `A` **including any
  explicit zeros**, matching Python `diffcp`. Covered by
  `tests/testthat/test-explicit-zeros.R`, including a finite-difference
  check that the recovered derivative is correct and not merely present.

## diffcp 0.1.1

CRAN release: 2026-05-27

- Fix `scs` issue with adaptive scaling triggered on R with *no long
  double* and also with MKL tolerance.
- Fix typo issues pointed out by CRAN.

## diffcp 0.1.0

CRAN release: 2026-05-19

First public release of the R port of
[`diffcp`](https://github.com/cvxgrp/diffcp).

### Features

- `solve_only(A, b, c, cone_dict, ...)` — forward solve of a convex cone
  program via Clarabel (default) or SCS, supporting all six standard
  cones (zero, nonneg, SOC, PSD, exponential, exponential dual).
  Quadratic objectives via the `P` argument are supported in
  forward-only mode (mirroring Python diffcp).

- `solve_and_derivative(A, b, c, cone_dict, mode, ...)` — solves the
  cone program and returns the optimal `(x, y, s)` together with two
  callables `D(dA, db, dc)` and `DT(dx, dy, ds)` that apply the
  derivative and its adjoint at the supplied perturbations.

  - `mode = "lsqr"` (default): matrix-free LSQR via the `M` operator.
  - `mode = "dense"`: dense `M` matrix factored with Eigen LDLT.

- `pi(x, cones, dual)` — projection onto a Cartesian product of cones,
  with `dual = TRUE` selecting the dual product cone.

- PSD-cone forward solves through Clarabel automatically permute the
  rows of `A` and entries of `b` from SCS lower-triangular ordering to
  Clarabel upper-triangular ordering, and permute `y` / `s` back on
  return.

### Implementation

- The numerical core is a faithful port of the C++ source in
  [cvxgrp/diffcp](https://github.com/cvxgrp/diffcp/tree/6a78143/cpp)
  (`linop.cpp`, `cones.cpp`, `deriv.cpp`, `lsqr.cpp`), called from R via
  `RcppEigen`.

- Each R file is annotated with `## DIFFCP SOURCE: <upstream path>`
  pointing at the corresponding lines in the upstream Python source.

### Tests

- 676 tests, faithful ports of every in-scope test in
  `tests/test_clarabel.py`, `tests/test_clarabel_psd.py`,
  `tests/test_scs.py`, and `tests/test_cone_prog_diff.py` from the
  upstream Python repository. Test fixtures are pinned to Python
  diffcp’s outputs at the same `atol` bounds Python uses (`1e-8` on
  finite-difference D-vs-FD agreement; `1e-4` on SCS-vs-Clarabel
  cross-checks; `1e-12` on round-trips).

### Out of scope (this release)

- The `lpgd` / `lpgd_left` / `lpgd_right` derivative modes (used by
  cvxpylayers for QP and degenerate-derivative cases).
- Batch parallel APIs (`solve_and_derivative_batch`,
  `solve_only_batch`).
- The `lsmr` mode.
- The ECOS solver branch (ECOS is unmaintained).
