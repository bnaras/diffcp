## CVXPY SOURCE: cvxpy/__init__.py (analogue: diffcp/__init__.py)
##
## Package-level documentation and namespace setup.

#' diffcp: Differentiating Through Cone Programs
#'
#' R port of the Python `diffcp` package (Agrawal, Barratt, Boyd 2019).
#' Computes the derivative of the optimal solution map of a convex cone
#' program with respect to its data, via the implicit-function theorem
#' on the homogeneous self-dual embedding.
#'
#' @name diffcp-package
#' @aliases diffcp
#' @keywords internal
#' @import Matrix
#' @useDynLib diffcp, .registration = TRUE
"_PACKAGE"

## CVXPY SOURCE: diffcp/cones.py:5-11
## Cone tag constants (string keys matching SCS / Clarabel conventions).

#' @export
ZERO     <- "z"
#' @export
EQ_DIM   <- ZERO
#' @export
POS      <- "l"
#' @export
SOC      <- "q"
#' @export
PSD      <- "s"
#' @export
EXP      <- "ep"
#' @export
EXP_DUAL <- "ed"

## The ordering of CONES matches SCS.
#' @export
CONES <- c(ZERO, POS, SOC, PSD, EXP, EXP_DUAL)
