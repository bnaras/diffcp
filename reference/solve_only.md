# Solve a cone program (forward only)

Mirrors `diffcp.cone_program.solve_only`. Solves minimize c^T x (+ 0.5
x^T P x) subject to A x + s = b, s in K and returns the optimal
`(x, y, s)`. Unlike `solve_and_derivative` this function does not build
the derivative closures.

## Usage

``` r
solve_only(
  A,
  b,
  c,
  cone_dict,
  warm_start = NULL,
  solve_method = "Clarabel",
  P = NULL,
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

- warm_start:

  Optional warm-start `list(x, y, s)`.

- solve_method:

  One of `"Clarabel"` (default) or `"SCS"`.

- P:

  Optional sparse `dgCMatrix` for QP objective.

- ...:

  Additional control parameters forwarded to the solver.

## Value

A list with elements `x`, `y`, `s`, `info`.
