## DIFFCP SOURCE: diffcp/__init__.py
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

## DIFFCP SOURCE: diffcp/cones.py:5-11

#' Cone tag constants
#'
#' String keys matching the SCS / Clarabel cone conventions. Used as
#' names in the `cone_dict` argument to [solve_only()] and
#' [solve_and_derivative()], and in the `(name, size)` pairs returned
#' by [parse_cone_dict()].
#'
#' @format Character scalars (and one character vector, `CONES`).
#' @details
#'   * `ZERO` (`"z"`): zero cone (equality constraints).
#'   * `EQ_DIM`: alias of `ZERO`.
#'   * `POS` (`"l"`): nonnegative orthant.
#'   * `SOC` (`"q"`): second-order cone.
#'   * `PSD` (`"s"`): symmetric positive semidefinite cone.
#'   * `EXP` (`"ep"`): primal exponential cone.
#'   * `EXP_DUAL` (`"ed"`): dual exponential cone.
#'   * `CONES`: the canonical SCS ordering of cone tags.
#' @name diffcp-cones
#' @keywords internal
NULL

#' @export
#' @rdname diffcp-cones
ZERO     <- "z"
#' @export
#' @rdname diffcp-cones
EQ_DIM   <- ZERO
#' @export
#' @rdname diffcp-cones
POS      <- "l"
#' @export
#' @rdname diffcp-cones
SOC      <- "q"
#' @export
#' @rdname diffcp-cones
PSD      <- "s"
#' @export
#' @rdname diffcp-cones
EXP      <- "ep"
#' @export
#' @rdname diffcp-cones
EXP_DUAL <- "ed"

#' @export
#' @rdname diffcp-cones
CONES <- c(ZERO, POS, SOC, PSD, EXP, EXP_DUAL)
