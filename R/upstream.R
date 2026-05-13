## R-SPECIFIC: upstream pin metadata helpers.
##
## The R port tracks an exact upstream Python `diffcp` commit. The
## canonical pin lives in `inst/UPSTREAM.dcf` (DCF so it's readable
## from R via `read.dcf`). Per-file `## DIFFCP SOURCE:` headers point
## at paths inside that pinned upstream tree.

#' Upstream Python diffcp pin
#'
#' Returns the upstream pin metadata for this R port: which
#' `cvxgrp/diffcp` version / commit the R sources track, the date of
#' that commit, and diffcp's own CVXPY runtime constraint. The
#' authoritative record lives in `inst/UPSTREAM.dcf` and is read at
#' call time.
#'
#' @return A named character vector with one element per DCF field
#'   (`Upstream`, `Version`, `Commit`, `ShortCommit`, `Date`, `URL`,
#'   `CVXPY`, `Snapshot`, `Notes`).
#' @export
#' @examples
#' upstream_info()
#' upstream_info()[["Version"]]
upstream_info <- function() {
  path <- system.file("UPSTREAM.dcf", package = "diffcp")
  if (!nzchar(path)) {
    stop("UPSTREAM.dcf not found in installed diffcp package.")
  }
  fields <- read.dcf(path)
  if (nrow(fields) != 1L) {
    stop("UPSTREAM.dcf must contain exactly one record.")
  }
  drop(fields)
}
