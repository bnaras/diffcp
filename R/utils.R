## CVXPY SOURCE: diffcp/cones.py (auxiliary utilities).
##
## Small helpers shared across cone projection / Jacobian code.

#' Symmetric vectorisation (SCS convention)
#'
#' Returns a vectorized representation of a symmetric matrix `X`,
#' with off-diagonal entries scaled by sqrt(2) to make the SCS-style
#' dot product `<vec_symm(A), vec_symm(B)> = trace(A B)` hold.
#'
#' Mirrors `diffcp.cones.vec_symm` (Python).
#'
#' @param X A symmetric matrix.
#' @returns A numeric vector of length `n*(n+1)/2`.
#' @export
vec_symm <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (ncol(X) != n) {
    cli::cli_abort("{.fn vec_symm} requires a square matrix; got {nrow(X)}x{ncol(X)}.")
  }
  Y <- X * sqrt(2)
  diag(Y) <- diag(X)
  ## Upper-triangular indices in row-major order (matches SCS).
  ## Python uses np.triu_indices(n) which returns (rows, cols) such that
  ## row <= col -- iterate column-major over upper triangle.
  idx <- which(upper.tri(Y, diag = TRUE), arr.ind = TRUE)
  ## Order so that we walk columns first within each row (i.e., for each
  ## row, all upper-triangle entries in that row in increasing col).
  ord <- order(idx[, 1L], idx[, 2L])
  Y[idx[ord, , drop = FALSE]]
}

#' Inverse of `vec_symm`
#'
#' @param x A numeric vector of length `n*(n+1)/2`.
#' @param dim The matrix dimension `n`.
#' @returns The corresponding `n x n` symmetric matrix.
#' @export
unvec_symm <- function(x, dim) {
  if (length(x) != dim * (dim + 1L) / 2L) {
    cli::cli_abort(
      "{.arg x} has length {length(x)}; expected {dim * (dim + 1L) / 2L}."
    )
  }
  X <- matrix(0, dim, dim)
  ## Walk upper triangle in row-major (matching vec_symm).
  k <- 1L
  for (i in seq_len(dim)) {
    for (j in i:dim) {
      X[i, j] <- x[k]
      k <- k + 1L
    }
  }
  X <- X + t(X)
  X[upper.tri(X)] <- X[upper.tri(X)] / sqrt(2)
  X[lower.tri(X)] <- X[lower.tri(X)] / sqrt(2)
  ## Diagonal already counted twice (once via X, once via t(X)); halve.
  diag(X) <- diag(X) / 2
  X
}

#' PSD cone vectorised dimension
#'
#' For a `dim x dim` symmetric matrix, returns the length of its
#' upper-triangular packed representation `dim*(dim+1)/2`.
#' @param dim Matrix dimension.
#' @returns Integer.
#' @noRd
.vec_psd_dim <- function(dim) {
  as.integer(dim * (dim + 1L) / 2L)
}

#' Inverse: from packed length back to matrix dimension.
#' @noRd
.psd_dim <- function(size) {
  ## solve k*(k+1)/2 = size -> k = (-1 + sqrt(1 + 8*size)) / 2
  k <- (-1 + sqrt(1 + 8 * size)) / 2
  as.integer(round(k))
}
