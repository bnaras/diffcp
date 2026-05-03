## DIFFCP SOURCE: diffcp/cpp/src/deriv.cpp + diffcp/cone_program.py:632-902
##
## Dense derivative path for `solve_and_derivative(mode = "dense")`.
##
## Mathematical setup (Busseti, Moursi, Boyd 2019):
##   z = (u, v, w) with u = x, v = y - s, w = 1.
##   Q = [[ 0,    A^T, c ];
##        [ -A,  0,    b ];
##        [ -c^T, -b^T, 0 ]]
##   pi(z) = (u, Pi_{K^*}(v), max(w, 0))   (skewed projection)
##   D pi(z) is block-diagonal with blocks (I_n, dPi_{K^*}(v),
##     1{w >= 0}); a (N x N) matrix where N = m + n + 1.
##   M = (Q - I) D pi(z) + I
##   The forward derivative solves M dz = dQ pi(z); the adjoint solves
##   M^T r = dz.  The Python source uses normal equations
##     (M^T M) dz = M^T rhs  (and (M M^T) r = M dz for the adjoint),
##   which gracefully handles boundary cases where M is singular.

# -- Skewed pi onto R^n x K^* x R_+ -------------------------------
## DIFFCP SOURCE: diffcp/cone_program.py:85-92 (the *other* `pi`).
.skew_pi <- function(z, cones, n_x, n_y) {
  u <- z[seq_len(n_x)]
  v <- z[(n_x + 1L):(n_x + n_y)]
  w <- z[n_x + n_y + 1L]
  c(u, pi(v, cones, dual = TRUE), max(w, 0))
}

# -- Build the (N x N) Q block ------------------------------------
## DIFFCP SOURCE: diffcp/cone_program.py:676-680 + 716-720 (used both
## for Q at the optimal z and for dQ in the derivative closure).
.build_Q <- function(A, b, c) {
  m <- nrow(A); n <- ncol(A)
  N <- m + n + 1L
  ## Q is dense for now; M_dense ultimately needs it dense anyway.
  Q <- matrix(0, N, N)
  AT <- as.matrix(Matrix::t(A))
  Aden <- as.matrix(A)
  Q[1L:n,                  (n + 1L):(n + m)]     <- AT
  Q[1L:n,                  N]                    <- c
  Q[(n + 1L):(n + m),      1L:n]                 <- -Aden
  Q[(n + 1L):(n + m),      N]                    <- b
  Q[N,                     1L:n]                 <- -c
  Q[N,                     (n + 1L):(n + m)]     <- -b
  Q
}

# -- D pi(z), dense ------------------------------------------------
## DIFFCP SOURCE: diffcp/cpp/src/deriv.cpp::dpi_dense (lines 30-42).
.dpi_dense <- function(u, v, w, cones) {
  n <- length(u); m <- length(v)
  N <- n + m + 1L
  D <- matrix(0, N, N)
  diag(D)[1L:n] <- 1
  ## Materialise the dual-cone projection Jacobian as an explicit (m,m)
  ## block by applying its matvec to each standard basis vector. For
  ## the cones we currently support (Zero, Nonneg, SOC) this is cheap
  ## and exact; PSD/EXP will need a denser path.
  Dproj <- .dprojection(v, cones, dual = TRUE)
  Mblock <- matrix(0, m, m)
  for (j in seq_len(m)) {
    e <- numeric(m); e[j] <- 1
    Mblock[, j] <- Dproj$matvec(e)
  }
  D[(n + 1L):(n + m), (n + 1L):(n + m)] <- Mblock
  D[N, N] <- if (w >= 0) 1 else 0
  D
}

# -- M dense -------------------------------------------------------
## DIFFCP SOURCE: diffcp/cpp/src/deriv.cpp::M_dense (lines 44-51).
.M_dense <- function(Q_dense, cones, u, v, w) {
  N <- length(u) + length(v) + 1L
  eyeN <- diag(N)
  (Q_dense - eyeN) %*% .dpi_dense(u, v, w, cones) + eyeN
}

# -- Normal-equation solves ---------------------------------------
## DIFFCP SOURCE: diffcp/cpp/src/deriv.cpp::_solve_derivative_dense
##                   (lines 53-57; M^T M dz = M^T rhs)
##                + diffcp/cpp/src/deriv.cpp::_solve_adjoint_derivative_dense
##                   (lines 59-63; M M^T r = M dz).
##
## The Python source caches no factorisation (see TODO in deriv.cpp).
##
## R-specific note: M is typically rank-deficient at active-set
## boundaries (e.g. an LP at optimality with most slacks zero), and
## both base R's `solve()` and `qr.solve()` fail in that case.
## Eigen's `.ldlt().solve(MT * rhs)` succeeds because LDLT on a
## positive-semidefinite matrix returns a particular solution rather
## than erroring. The closest general-purpose R analogue is the
## Moore-Penrose pseudoinverse via SVD (`.pinv_solve` below), which
## both handles rank deficiency and gives the minimum-norm solution.
.pinv_solve <- function(A, b, tol = NULL) {
  s <- svd(A)
  if (is.null(tol)) {
    tol <- max(dim(A)) * .Machine$double.eps * max(s$d)
  }
  d_inv <- ifelse(s$d > tol, 1 / s$d, 0)
  s$v %*% (d_inv * crossprod(s$u, b))
}
.solve_derivative_dense <- function(M, MT, rhs) {
  .pinv_solve(M, rhs)
}
.solve_adjoint_derivative_dense <- function(M, MT, dz) {
  .pinv_solve(MT, dz)
}

# -- Build dQ from (dA, db, dc) -----------------------------------
## DIFFCP SOURCE: diffcp/cone_program.py:716-720 (inside derivative()).
.build_dQ <- function(dA, db, dc) {
  m <- nrow(dA); n <- ncol(dA)
  N <- m + n + 1L
  dQ <- matrix(0, N, N)
  dAT <- as.matrix(Matrix::t(dA))
  dAden <- as.matrix(dA)
  dQ[1L:n,                  (n + 1L):(n + m)] <- dAT
  dQ[1L:n,                  N]                <- dc
  dQ[(n + 1L):(n + m),      1L:n]             <- -dAden
  dQ[(n + 1L):(n + m),      N]                <- db
  dQ[N,                     1L:n]             <- -dc
  dQ[N,                     (n + 1L):(n + m)] <- -db
  dQ
}

# -- Build the closures used by `solve_and_derivative` ------------
## DIFFCP SOURCE: diffcp/cone_program.py:696-774 (derivative +
## adjoint_derivative closures).
##
## NOTE on signatures: Python returns `dA` as a SciPy CSC sparse
## matrix with the sparsity pattern of the original `A`. R's analogue
## here is a `dgCMatrix` built from the same `(rows, cols)`.
.make_derivative_closures <- function(A, b, c, x, y, s, cones, mode) {
  if (!mode %in% c("dense", "lsqr")) {
    cli::cli_abort("Unsupported mode {.val {mode}}; supported: 'dense', 'lsqr'.")
  }
  m <- nrow(A); n <- ncol(A); N <- m + n + 1L

  ## Capture nonzero pattern of A *before* we drop zeros, matching
  ## diffcp/cone_program.py:643-653 (NaN-trick to record true zeros).
  A_drop <- Matrix::drop0(A)
  A_csc  <- methods::as(A_drop, "CsparseMatrix")
  nz <- Matrix::summary(A_csc)
  rows <- nz$i  # 1-based
  cols <- nz$j  # 1-based

  z <- c(x, y - s, 1)
  Q <- .build_Q(A_csc, b, c)
  pi_z <- .skew_pi(z, cones, n_x = n, n_y = m)

  ## D_proj_dual_cone: Jacobian of Pi_{K^*} at v = y - s.
  Dproj_dual <- .dprojection(z[(n + 1L):(n + m)], cones, dual = TRUE)

  if (mode == "dense") {
    M  <- .M_dense(Q, cones, u = x, v = y - s, w = 1)
    MT <- t(M)
  } else {
    cli::cli_abort(c(
      "Mode {.val lsqr} not yet implemented in this development version.",
      "i" = "Use {.code mode = \"dense\"} for now; LSQR lands in a follow-on commit."
    ))
  }

  derivative <- function(dA, db, dc) {
    dA_csc <- methods::as(methods::as(dA, "CsparseMatrix"), "dMatrix")
    dQ  <- .build_dQ(dA_csc, db, dc)
    rhs <- as.numeric(dQ %*% pi_z)
    if (isTRUE(all.equal(rhs, numeric(N), tolerance = 0)) ||
        all(abs(rhs) < .Machine$double.eps)) {
      dz <- numeric(N)
    } else {
      dz <- as.numeric(.solve_derivative_dense(M, MT, rhs))
    }
    du <- dz[1L:n]
    dv <- dz[(n + 1L):(n + m)]
    dw <- dz[N]
    dx <- du - x * dw
    dy <- Dproj_dual$matvec(dv) - y * dw
    ds <- Dproj_dual$matvec(dv) - dv - s * dw
    list(dx = -dx, dy = -dy, ds = -ds)
  }

  adjoint_derivative <- function(dx, dy, ds) {
    dw <- -(sum(x * dx) + sum(y * dy) + sum(s * ds))
    dz <- c(dx, Dproj_dual$rmatvec(dy + ds) - ds, dw)
    if (all(abs(dz) < .Machine$double.eps)) {
      r <- numeric(N)
    } else {
      r <- as.numeric(.solve_adjoint_derivative_dense(M, MT, dz))
    }
    ## Python uses 0-based (rows, cols); we use 1-based here so the
    ## indices `rows + n` become `rows + n` directly in 1-based terms.
    values <- pi_z[cols] * r[rows + n] - pi_z[n + rows] * r[cols]
    dA <- Matrix::sparseMatrix(i = rows, j = cols, x = values,
                               dims = c(m, n))
    db <- pi_z[(n + 1L):(n + m)] * r[N] - pi_z[N] * r[(n + 1L):(n + m)]
    dc <- pi_z[1L:n]              * r[N] - pi_z[N] * r[1L:n]
    list(dA = dA, db = db, dc = dc)
  }

  list(D = derivative, DT = adjoint_derivative)
}
