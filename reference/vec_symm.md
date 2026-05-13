# Symmetric vectorization (SCS convention)

Returns a vectorized representation of a symmetric matrix `X`, with
off-diagonal entries scaled by sqrt(2) to make the SCS-style dot product
`<vec_symm(A), vec_symm(B)> = trace(A B)` hold.

## Usage

``` r
vec_symm(X)
```

## Arguments

- X:

  A symmetric matrix.

## Value

A numeric vector of length `n*(n+1)/2`.

## Details

Mirrors `diffcp.cones.vec_symm` (Python).
