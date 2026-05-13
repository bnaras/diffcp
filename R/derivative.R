## DIFFCP SOURCE: diffcp/cone_program.py:632-902 @ v1.1.8 (6a78143)
## Upstream pin: see inst/UPSTREAM.dcf / upstream_info().
##
## R-side orchestration of `solve_and_derivative`. All numerical heavy
## lifting (cone projection Jacobians, M operator, dense LDLT solve,
## matrix-free LSQR) is delegated to C++ via the wrappers in
## src/wrapper.cpp.
##
## This file mirrors the Python source line-for-line, except that the
## per-iteration linear-algebra routines (`dpi_dense`, `M_dense`,
## `_solve_derivative_dense`, `_solve_adjoint_derivative_dense`,
## `M_operator`, `lsqr`) are C++ functions reachable from R.

# -- Skewed pi onto R^n x K^* x R_+ -------------------------------
## DIFFCP SOURCE: diffcp/cone_program.py:85-92 (`pi(z, cones)`).
.skew_pi <- function(z, cones, n_x, n_y) {
  u <- z[seq_len(n_x)]
  v <- z[(n_x + 1L):(n_x + n_y)]
  w <- z[n_x + n_y + 1L]
  c(u, pi(v, cones, dual = TRUE), max(w, 0))
}

# -- Build the (N x N) Q block ------------------------------------
## DIFFCP SOURCE: diffcp/cone_program.py:676-680 (Q at z*) and
## :716-720 (dQ inside derivative()). Returns a sparse dgCMatrix.
.build_Q <- function(A, b, c) {
  m <- nrow(A); n <- ncol(A)
  N <- m + n + 1L
  ## Sparse-bmat equivalent: explicitly build (i, j, x) triplets.
  AT <- methods::as(Matrix::t(A), "TsparseMatrix")
  Acoo <- methods::as(A, "TsparseMatrix")
  ## block (1:n, (n+1):(n+m)) = A^T  (rows of A^T are A's columns)
  i_at <- AT@i + 1L
  j_at <- AT@j + 1L + n
  ## block ((n+1):(n+m), 1:n) = -A
  i_a  <- Acoo@i + 1L + n
  j_a  <- Acoo@j + 1L
  ## column N: c stacked over b stacked over 0; row N: -c^T stacked
  ## with -b^T stacked with 0 (mirroring deriv.cpp / cone_program.py).
  i_col <- c(seq_len(n), seq.int(n + 1L, n + m))
  j_col <- rep.int(N, n + m)
  x_col <- c(c, b)
  i_row <- rep.int(N, n + m)
  j_row <- c(seq_len(n), seq.int(n + 1L, n + m))
  x_row <- c(-c, -b)

  Matrix::sparseMatrix(
    i = c(i_at, i_a, i_col, i_row),
    j = c(j_at, j_a, j_col, j_row),
    x = c(AT@x, -Acoo@x, x_col, x_row),
    dims = c(N, N)
  )
}

# -- Build the closures used by `solve_and_derivative` ------------
## DIFFCP SOURCE: diffcp/cone_program.py:696-774 (derivative +
## adjoint_derivative closures).
.make_derivative_closures <- function(A, b, c, x, y, s, cones, mode) {
  if (!mode %in% c("dense", "lsqr")) {
    cli::cli_abort("Unsupported mode {.val {mode}}; supported: 'dense', 'lsqr'.")
  }
  m <- nrow(A); n <- ncol(A); N <- m + n + 1L

  A_drop <- Matrix::drop0(A)
  A_csc  <- methods::as(A_drop, "CsparseMatrix")
  nz <- Matrix::summary(A_csc)
  rows <- nz$i  # 1-based
  cols <- nz$j  # 1-based

  z <- c(x, y - s, 1)
  Q_sp <- .build_Q(A_csc, b, c)
  pi_z <- .skew_pi(z, cones, n_x = n, n_y = m)

  ## D pi(v) for the dual cone, materialised once via the C++ side so
  ## the same code path serves all cones (Zero/Nonneg/SOC today,
  ## PSD/EXP when those land).
  Dproj_dual <- cpp_dprojection_dense(z[(n + 1L):(n + m)], cones, dual = TRUE)

  if (mode == "dense") {
    Q_dense <- as.matrix(Q_sp)
    M  <- cpp_M_dense(Q_dense, cones, u = x, v = y - s, w = 1)
    MT <- t(M)
  }

  derivative <- function(dA, db, dc) {
    dA_csc <- methods::as(dA, "CsparseMatrix")
    dQ_sp  <- .build_Q(dA_csc, db, dc)
    rhs    <- as.numeric(dQ_sp %*% pi_z)
    if (all(abs(rhs) < .Machine$double.eps)) {
      dz <- numeric(N)
    } else if (mode == "dense") {
      dz <- as.numeric(cpp_solve_derivative_dense(M, MT, rhs))
    } else {
      Q_map <- methods::as(Q_sp, "CsparseMatrix")
      out <- cpp_lsqr_M(Q_map, cones, x, y - s, 1, rhs, transpose = FALSE)
      dz <- as.numeric(out$solution)
    }
    du <- dz[1L:n]
    dv <- dz[(n + 1L):(n + m)]
    dw <- dz[N]
    Ddv <- as.numeric(Dproj_dual %*% dv)
    dx <- du - x * dw
    dy <- Ddv - y * dw
    ds <- Ddv - dv - s * dw
    list(dx = -dx, dy = -dy, ds = -ds)
  }

  adjoint_derivative <- function(dx, dy, ds) {
    dw <- -(sum(x * dx) + sum(y * dy) + sum(s * ds))
    dz <- c(dx, as.numeric(crossprod(Dproj_dual, dy + ds)) - ds, dw)
    if (all(abs(dz) < .Machine$double.eps)) {
      r <- numeric(N)
    } else if (mode == "dense") {
      r <- as.numeric(cpp_solve_adjoint_derivative_dense(M, MT, dz))
    } else {
      Q_map <- methods::as(Q_sp, "CsparseMatrix")
      out <- cpp_lsqr_M(Q_map, cones, x, y - s, 1, dz, transpose = TRUE)
      r <- as.numeric(out$solution)
    }
    values <- pi_z[cols] * r[rows + n] - pi_z[n + rows] * r[cols]
    dA <- Matrix::sparseMatrix(i = rows, j = cols, x = values,
                               dims = c(m, n))
    db <- pi_z[(n + 1L):(n + m)] * r[N] - pi_z[N] * r[(n + 1L):(n + m)]
    dc <- pi_z[1L:n]              * r[N] - pi_z[N] * r[1L:n]
    list(dA = dA, db = db, dc = dc)
  }

  list(D = derivative, DT = adjoint_derivative)
}
