# Projection onto a Cartesian product of cones

Mirrors `diffcp.cones.pi`.

## Usage

``` r
pi(x, cones, dual = FALSE)
```

## Arguments

- x:

  A numeric vector.

- cones:

  A list of `list(name, size)` pairs, as produced by `parse_cone_dict`.

- dual:

  If `TRUE`, project onto the dual cones.

## Value

A numeric vector of the same length as `x`.
