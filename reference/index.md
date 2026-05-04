# Package index

## Solve a cone program

Forward solve and the differentiable solve. Both accept (A, b, c) in SCS
standard form and a `cone_dict` of cone sizes.

- [`solve_only()`](https://bnaras.github.io/diffcp/reference/solve_only.md)
  : Solve a cone program (forward only)
- [`solve_and_derivative()`](https://bnaras.github.io/diffcp/reference/solve_and_derivative.md)
  : Solve a cone program and return forward / adjoint derivative
  operators

## Cone projection

Projection onto a Cartesian product of cones, plus the SCS-convention
vec / unvec helpers for symmetric matrices.

- [`pi()`](https://bnaras.github.io/diffcp/reference/pi.md) : Projection
  onto a Cartesian product of cones

- [`vec_symm()`](https://bnaras.github.io/diffcp/reference/vec_symm.md)
  : Symmetric vectorisation (SCS convention)

- [`unvec_symm()`](https://bnaras.github.io/diffcp/reference/unvec_symm.md)
  :

  Inverse of `vec_symm`

- [`parse_cone_dict()`](https://bnaras.github.io/diffcp/reference/parse_cone_dict.md)
  : Parse an SCS-style cone dictionary into an ordered list

## Cone tag constants

- [`ZERO`](https://bnaras.github.io/diffcp/reference/diffcp-cones.md)
  [`EQ_DIM`](https://bnaras.github.io/diffcp/reference/diffcp-cones.md)
  [`POS`](https://bnaras.github.io/diffcp/reference/diffcp-cones.md)
  [`SOC`](https://bnaras.github.io/diffcp/reference/diffcp-cones.md)
  [`PSD`](https://bnaras.github.io/diffcp/reference/diffcp-cones.md)
  [`EXP`](https://bnaras.github.io/diffcp/reference/diffcp-cones.md)
  [`EXP_DUAL`](https://bnaras.github.io/diffcp/reference/diffcp-cones.md)
  [`CONES`](https://bnaras.github.io/diffcp/reference/diffcp-cones.md) :
  Cone tag constants

## Package

- [`diffcp-package`](https://bnaras.github.io/diffcp/reference/diffcp-package.md)
  [`diffcp`](https://bnaras.github.io/diffcp/reference/diffcp-package.md)
  : diffcp: Differentiating Through Cone Programs
