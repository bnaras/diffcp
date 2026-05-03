# Cone tag constants

String keys matching the SCS / Clarabel cone conventions. Used as names
in the `cone_dict` argument to
[`solve_only()`](https://bnaras.github.io/diffcp/reference/solve_only.md)
and
[`solve_and_derivative()`](https://bnaras.github.io/diffcp/reference/solve_and_derivative.md),
and in the `(name, size)` pairs returned by
[`parse_cone_dict()`](https://bnaras.github.io/diffcp/reference/parse_cone_dict.md).

## Usage

``` r
ZERO

EQ_DIM

POS

SOC

PSD

EXP

EXP_DUAL

CONES
```

## Format

Character scalars (and one character vector, `CONES`).

An object of class `character` of length 1.

An object of class `character` of length 1.

An object of class `character` of length 1.

An object of class `character` of length 1.

An object of class `character` of length 1.

An object of class `character` of length 1.

An object of class `character` of length 1.

An object of class `character` of length 6.

## Details

- `ZERO` (`"z"`): zero cone (equality constraints).

- `EQ_DIM`: alias of `ZERO`.

- `POS` (`"l"`): nonnegative orthant.

- `SOC` (`"q"`): second-order cone.

- `PSD` (`"s"`): symmetric positive semidefinite cone.

- `EXP` (`"ep"`): primal exponential cone.

- `EXP_DUAL` (`"ed"`): dual exponential cone.

- `CONES`: the canonical SCS ordering of cone tags.
