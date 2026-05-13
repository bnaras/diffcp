# Upstream Python diffcp pin

Returns the upstream pin metadata for this R port: which `cvxgrp/diffcp`
version / commit the R sources track, the date of that commit, and
diffcp's own CVXPY runtime constraint. The authoritative record lives in
`inst/UPSTREAM.dcf` and is read at call time.

## Usage

``` r
upstream_info()
```

## Value

A named character vector with one element per DCF field (`Upstream`,
`Version`, `Commit`, `ShortCommit`, `Date`, `URL`, `CVXPY`, `Snapshot`,
`Notes`).

## Examples

``` r
upstream_info()
#>                                                                                                                                                                                                                                                                                                                                                                                 Upstream 
#>                                                                                                                                                                                                                                                                                                                                                                          "cvxgrp/diffcp" 
#>                                                                                                                                                                                                                                                                                                                                                                                  Version 
#>                                                                                                                                                                                                                                                                                                                                                                                  "1.1.8" 
#>                                                                                                                                                                                                                                                                                                                                                                                   Commit 
#>                                                                                                                                                                                                                                                                                                                                               "6a78143daccc80db4c27da39ba7edf5c9a388586" 
#>                                                                                                                                                                                                                                                                                                                                                                              ShortCommit 
#>                                                                                                                                                                                                                                                                                                                                                                                "6a78143" 
#>                                                                                                                                                                                                                                                                                                                                                                                     Date 
#>                                                                                                                                                                                                                                                                                                                                                                             "2026-01-26" 
#>                                                                                                                                                                                                                                                                                                                                                                                      URL 
#>                                                                                                                                                                                                                                                                                                                                                       "https://github.com/cvxgrp/diffcp" 
#>                                                                                                                                                                                                                                                                                                                                                                                    CVXPY 
#>                                                                                                                                                                                                                                                                                                                                                                               ">= 1.6.3" 
#>                                                                                                                                                                                                                                                                                                                                                                                 Snapshot 
#>                                                                                                                                                                                                                                                                                                                                   "/Users/naras/research/cvxr/new_design/diffing/diffcp" 
#>                                                                                                                                                                                                                                                                                                                                                                                    Notes 
#> "Pinned snapshot of the Python diffcp source that this R port\ntracks. Per-file `## DIFFCP SOURCE:` headers reference paths inside\nthis upstream tree. The CVXPY constraint is diffcp's own runtime\ndependency on the Python side; the R port has no direct CVXPY\ndependency but ships into a CVXR ecosystem aligned to a CVXPY\nversion (see CVXR's own DESCRIPTION/Config fields)." 
upstream_info()[["Version"]]
#> [1] "1.1.8"
```
