# Parse an SCS-style cone dictionary into an ordered list

Mirrors `diffcp.cones.parse_cone_dict`.

## Usage

``` r
parse_cone_dict(cone_dict)
```

## Arguments

- cone_dict:

  A named list with keys among `"z"`, `"l"`, `"q"`, `"s"`, `"ep"`,
  `"ed"`. Values are scalar dimensions (for `z`, `l`, `ep`, `ed`) or
  integer vectors of dimensions (for `q`, `s`).

## Value

A list of `list(name, size)` pairs, in the canonical SCS cone order
(`CONES`).
