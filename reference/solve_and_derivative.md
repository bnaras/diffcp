# Solve a cone program and return forward / adjoint derivative operators

Mirrors `diffcp.cone_program.solve_and_derivative`. Solves minimize c^T
x s.t. A x + s = b, s in K (with optional QP `0.5 x^T P x` term) and
returns the optimal `(x, y, s)` together with closures `D` (forward) and
`DT` (adjoint) that map perturbations of `(A, b, c, [P])` to
perturbations of `(x, y, s)` and vice versa.

## Usage

``` r
solve_and_derivative(
  A,
  b,
  c,
  cone_dict,
  P = NULL,
  solve_method = "Clarabel",
  mode = "lsqr",
  warm_start = NULL,
  ...
)
```

## Arguments

- A:

  A sparse `dgCMatrix` constraint matrix.

- b:

  A numeric offset vector.

- c:

  A numeric objective coefficient vector.

- cone_dict:

  A named list with cone sizes (keys among `"z"`, `"l"`, `"q"`, `"s"`,
  `"ep"`, `"ed"`).

- P:

  Optional sparse `dgCMatrix` for QP objective.

- solve_method:

  One of `"Clarabel"` (default) or `"SCS"`.

- mode:

  Differentiation mode: `"lsqr"` (default), `"dense"`, or `"lpgd"`.

- warm_start:

  Optional warm-start `list(x, y, s)`.

- ...:

  Additional control parameters forwarded to the solver.

## Value

A list with elements `x`, `y`, `s`, `info`, `D`, `DT`. `info` is the
solver-status block returned by the underlying solver (`status`, `iter`,
`solveTime`, `pobj`, ...); see `solve_only` for the same shape.
