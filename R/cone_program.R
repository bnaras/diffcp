## DIFFCP SOURCE: diffcp/cone_program.py
##
## Forward solve dispatch (Clarabel / SCS) plus the public
## `solve_only` and `solve_and_derivative` entry points. Derivative
## computation (M operator, LSQR, dense modes) lands in subsequent
## commits; for now `solve_and_derivative` errors with a clear
## phase-pointer message, while `solve_only` is fully functional.

# -- PSD permutation: SCS lower-triangular <-> Clarabel upper-triangular --
##
## DIFFCP SOURCE: diffcp/cone_program.py:14-82 (`permute_psd_rows`,
## `inverse_permute_psd_solution`).
##
## SCS packs a vec(X) for an n x n PSD matrix as the lower triangle in
## column-major order: (X11, X21, ..., Xn1, X22, X32, ..., Xnn).
## Clarabel packs the upper triangle in column-major order:
## (X11, X12, X22, X13, X23, X33, ...).  For symmetric X both encode
## the same matrix, but the *positions* differ, so the rows of A and
## entries of b corresponding to one PSD block must be permuted before
## handing to Clarabel; y and s come back in Clarabel order and must
## be permuted back to SCS order.
##
## We compute two permutations matching Python:
##   post[k] = clarabel index for the entry whose SCS index is k
##   pre[k]  = scs index for the entry whose clarabel index is k
## (post and pre are inverses of each other).
.psd_perms <- function(n) {
  ## np.tril_indices(n): for r=0..n-1, for c=0..r yield (r, c)
  ## np.triu_indices(n): for r=0..n-1, for c=r..n-1 yield (r, c)
  tril <- do.call(rbind, lapply(0:(n - 1L),
                                function(r) cbind(r, 0:r)))
  triu <- do.call(rbind, lapply(0:(n - 1L),
                                function(r) cbind(r, r:(n - 1L))))
  ## ravel_multi_index((cols, rows), (n, n)) row-major:
  ##   linear = cols * n + rows
  tril_multi <- tril[, 2L] * n + tril[, 1L]
  triu_multi <- triu[, 2L] * n + triu[, 1L]
  ## np.argsort -> R order(); convert to 0-based offsets.
  list(
    post = order(tril_multi) - 1L,   # SCS k -> Clarabel position
    pre  = order(triu_multi) - 1L    # Clarabel k -> SCS position
  )
}

## Permute rows of a CsparseMatrix A and entries of b for ONE PSD
## block of dimension `n` (n*(n+1)/2 vectorised entries) starting at
## 1-based row index `row_offset`.
.permute_psd_rows <- function(A, b, n, row_offset) {
  perms <- .psd_perms(n)
  n_rows <- length(perms$post)

  At <- methods::as(A, "TsparseMatrix")
  rows <- At@i + 1L                    # 1-based
  cols <- At@j + 1L
  vals <- At@x
  in_block <- rows >= row_offset & rows < row_offset + n_rows
  rows[in_block] <- row_offset + perms$post[rows[in_block] - row_offset + 1L]

  new_A <- Matrix::sparseMatrix(i = rows, j = cols, x = vals, dims = dim(A))
  new_b <- b
  new_b[row_offset:(row_offset + n_rows - 1L)] <-
    b[row_offset + perms$pre + 1L - 1L]   # 1-based
  list(A = new_A, b = new_b)
}

## Inverse permutation applied to y and s after a Clarabel solve.
.inverse_permute_psd_solution <- function(y, s, n, row_offset) {
  perms <- .psd_perms(n)
  n_rows <- length(perms$pre)
  idx <- row_offset + perms$pre   # 1-based source positions
  new_y <- y
  new_s <- s
  new_y[row_offset:(row_offset + n_rows - 1L)] <- y[idx]
  new_s[row_offset:(row_offset + n_rows - 1L)] <- s[idx]
  list(y = new_y, s = new_s)
}

# -- Internal: dispatch the forward solve to a concrete solver ----
## DIFFCP SOURCE: diffcp/cone_program.py:397-630 (solve_internal).
##
## Differences from the Python source in this development version:
##   * ECOS branch is omitted (R `ECOSolveR` API differs; deferred
##     until we need ECOS-specific behaviour for tests).
##   * Clarabel PSD row-permutation (`permute_psd_rows` /
##     `inverse_permute_psd_solution`) is omitted — PSD support
##     arrives in Phase 2b alongside the PSD cone projection.
##   * SCS `acceleration_lookback=0` retry on inaccurate status is
##     dropped (R `scs` package takes the same control name and we
##     can add the retry once we hit an inaccurate case in tests).
.solve_internal <- function(A, b, c, cone_dict,
                            solve_method = NULL,
                            warm_start = NULL,
                            raise_on_error = TRUE,
                            P = NULL,
                            ...) {
  kwargs <- list(...)

  if (is.null(solve_method)) {
    psd_cone <- !is.null(cone_dict[["s"]]) && length(cone_dict[["s"]]) > 0L
    ed_cone  <- !is.null(cone_dict[["ed"]]) && cone_dict[["ed"]] > 0L
    solve_method <- if (psd_cone || ed_cone) "SCS" else "Clarabel"
  }

  if (solve_method == "Clarabel" || solve_method == "CLARABEL") {
    if (!requireNamespace("clarabel", quietly = TRUE)) {
      cli::cli_abort("The {.pkg clarabel} package is required for Clarabel.")
    }
    if (!is.null(warm_start)) {
      cli::cli_abort("Clarabel does not support warm starting.")
    }

    ## Clarabel uses the same SCS-style cones list. Translate `f`->`z`
    ## (free cone alias) and drop `ed` (Clarabel does not implement
    ## dual exponential cones natively).  PSD cones need an additional
    ## row permutation on A and b, applied below per cone in cone-order
    ## start_row = z + f + l + sum(q); then for each PSD of size v,
    ## permute and advance by v*(v+1)/2.
    cones_arg <- list()
    if (!is.null(cone_dict[["z"]]) && cone_dict[["z"]] > 0L)
      cones_arg$z <- as.integer(cone_dict[["z"]])
    if (!is.null(cone_dict[["f"]]) && cone_dict[["f"]] > 0L) {
      cones_arg$z <- (cones_arg$z %||% 0L) + as.integer(cone_dict[["f"]])
    }
    if (!is.null(cone_dict[["l"]]) && cone_dict[["l"]] > 0L)
      cones_arg$l <- as.integer(cone_dict[["l"]])
    if (!is.null(cone_dict[["q"]]) && length(cone_dict[["q"]]) > 0L)
      cones_arg$q <- as.integer(cone_dict[["q"]])
    if (!is.null(cone_dict[["s"]]) && length(cone_dict[["s"]]) > 0L)
      cones_arg$s <- as.integer(cone_dict[["s"]])
    if (!is.null(cone_dict[["ep"]]) && cone_dict[["ep"]] > 0L)
      cones_arg$ep <- as.integer(cone_dict[["ep"]])

    ## Apply PSD row permutations to A, b in SCS->Clarabel direction.
    if (!is.null(cone_dict[["s"]]) && length(cone_dict[["s"]]) > 0L) {
      start_row <- 1L  # 1-based
      if (!is.null(cone_dict[["z"]])) start_row <- start_row + as.integer(cone_dict[["z"]])
      if (!is.null(cone_dict[["f"]])) start_row <- start_row + as.integer(cone_dict[["f"]])
      if (!is.null(cone_dict[["l"]])) start_row <- start_row + as.integer(cone_dict[["l"]])
      if (!is.null(cone_dict[["q"]])) start_row <- start_row + sum(as.integer(cone_dict[["q"]]))
      for (v in as.integer(cone_dict[["s"]])) {
        permuted <- .permute_psd_rows(A, b, v, start_row)
        A <- permuted$A
        b <- permuted$b
        start_row <- start_row + (v * (v + 1L)) %/% 2L
      }
    }

    if (is.null(P)) {
      n <- length(c)
      P <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                dims = c(n, n))
    }

    control <- do.call(clarabel::clarabel_control, c(list(verbose = FALSE), kwargs))
    A_csc <- methods::as(A, "CsparseMatrix")
    P_csc <- methods::as(P, "CsparseMatrix")
    sol <- clarabel::clarabel(A = A_csc, b = b, q = c, P = P_csc,
                              cones = cones_arg, control = control)

    ## R clarabel: `sol$status` is an integer index into
    ## `solver_status_descriptions()`; the *names* of that vector are
    ## the canonical status tags ("Solved", "PrimalInfeasible", ...)
    ## that match Python clarabel's `str(solution.status)`.
    status_str  <- names(clarabel::solver_status_descriptions())[sol$status]
    CLARABEL2SCS_STATUS_MAP <- c(
      "Solved"                  = "Solved",
      "PrimalInfeasible"        = "Infeasible",
      "DualInfeasible"          = "Unbounded",
      "AlmostSolved"            = "Optimal Inaccurate",
      "AlmostPrimalInfeasible"  = "Infeasible Inaccurate",
      "AlmostDualInfeasible"    = "Unbounded Inaccurate"
    )
    scs_status <- CLARABEL2SCS_STATUS_MAP[status_str]
    if (is.na(scs_status)) scs_status <- "Failure"

    inaccurate_status <- c("Optimal Inaccurate", "Infeasible Inaccurate",
                           "Unbounded Inaccurate")
    if (scs_status == "Optimal Inaccurate") {
      ## Python's Clarabel branch (cone_program.py:540-626) does NOT
      ## error on AlmostSolved; downstream code accepts the answer.
      ## We only warn here.
      cli::cli_warn("Clarabel returned {.val {status_str}}.")
    } else if (scs_status != "Solved") {
      if (raise_on_error) {
        cli::cli_abort("Clarabel returned status {.val {status_str}}.")
      } else {
        return(list(x = sol$x, y = sol$z, s = sol$s,
                    info = list(status = scs_status), D = NULL, DT = NULL))
      }
    }

    y_out <- as.numeric(sol$z)
    s_out <- as.numeric(sol$s)

    ## Permute y and s back from Clarabel (upper-triangular) to SCS
    ## (lower-triangular) for each PSD cone block.
    if (!is.null(cone_dict[["s"]]) && length(cone_dict[["s"]]) > 0L) {
      start_row <- 1L
      if (!is.null(cone_dict[["z"]])) start_row <- start_row + as.integer(cone_dict[["z"]])
      if (!is.null(cone_dict[["f"]])) start_row <- start_row + as.integer(cone_dict[["f"]])
      if (!is.null(cone_dict[["l"]])) start_row <- start_row + as.integer(cone_dict[["l"]])
      if (!is.null(cone_dict[["q"]])) start_row <- start_row + sum(as.integer(cone_dict[["q"]]))
      for (v in as.integer(cone_dict[["s"]])) {
        permuted <- .inverse_permute_psd_solution(y_out, s_out, v, start_row)
        y_out <- permuted$y
        s_out <- permuted$s
        start_row <- start_row + (v * (v + 1L)) %/% 2L
      }
    }

    return(list(
      x = as.numeric(sol$x),
      y = y_out,
      s = s_out,
      info = list(
        status    = scs_status,
        solveTime = sol$solve_time %||% NA_real_,
        setupTime = -1,
        iter      = sol$iterations %||% NA_integer_,
        pobj      = sol$obj_val    %||% NA_real_
      )
    ))
  }

  if (solve_method == "SCS") {
    if (!requireNamespace("scs", quietly = TRUE)) {
      cli::cli_abort("The {.pkg scs} package is required for SCS.")
    }

    ## R `scs::scs` uses `eps_abs` / `eps_rel` (matching SCS >= 3.x).
    ## Mirror the Python convenience: a single `eps` becomes both.
    if (!is.null(kwargs[["eps"]])) {
      kwargs[["eps_abs"]] <- kwargs[["eps"]]
      kwargs[["eps_rel"]] <- kwargs[["eps"]]
      kwargs[["eps"]] <- NULL
    }
    kwargs$verbose <- kwargs$verbose %||% FALSE
    control <- do.call(scs::scs_control, kwargs)

    initial <- NULL
    if (!is.null(warm_start)) {
      initial <- list(x = warm_start[[1L]], y = warm_start[[2L]], s = warm_start[[3L]])
    }

    A_csc <- methods::as(A, "CsparseMatrix")
    ## R `scs` accepts P only as base `matrix` or
    ## `slam::simple_triplet_matrix` (no Matrix-package CsparseMatrix
    ## method).  Densify if sparse — fine in practice because P is
    ## typically small.
    P_in <- if (is.null(P)) NULL else as.matrix(P)
    sol <- scs::scs(A = A_csc, b = b, obj = c, P = P_in,
                    cone = cone_dict, initial = initial, control = control)

    status <- sol$info$status
    inaccurate_status <- c("Solved/Inaccurate",
                           "solved (inaccurate - reached max_iters)",
                           "solved (inaccurate - reached time_limit_secs)")
    if (status %in% inaccurate_status) {
      cli::cli_warn("Solved/Inaccurate.")
    } else if (tolower(status) != "solved") {
      if (raise_on_error) {
        cli::cli_abort("SCS returned status {.val {status}}.")
      } else {
        return(list(x = sol$x, y = sol$y, s = sol$s,
                    info = sol$info, D = NULL, DT = NULL))
      }
    }

    return(list(
      x = as.numeric(sol$x),
      y = as.numeric(sol$y),
      s = as.numeric(sol$s),
      info = sol$info
    ))
  }

  cli::cli_abort("Solver {.val {solve_method}} not supported.")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Solve a cone program (forward only)
#'
#' Mirrors `diffcp.cone_program.solve_only`. Solves
#'   minimise  c^T x  (+ 0.5 x^T P x)
#'   subject to  A x + s = b,   s in K
#' and returns the optimal `(x, y, s)`. Unlike `solve_and_derivative`
#' this function does not build the derivative closures.
#'
#' @param A A sparse `dgCMatrix` constraint matrix.
#' @param b A numeric offset vector.
#' @param c A numeric objective coefficient vector.
#' @param cone_dict A named list with cone sizes (keys among
#'   `"z"`, `"l"`, `"q"`, `"s"`, `"ep"`, `"ed"`).
#' @param warm_start Optional warm-start `list(x, y, s)`.
#' @param solve_method One of `"Clarabel"` (default) or `"SCS"`.
#' @param P Optional sparse `dgCMatrix` for QP objective.
#' @param ... Additional control parameters forwarded to the solver.
#' @returns A list with elements `x`, `y`, `s`, `info`.
#' @export
solve_only <- function(A, b, c, cone_dict,
                       warm_start = NULL,
                       solve_method = "Clarabel",
                       P = NULL,
                       ...) {
  if (any(is.nan(A@x))) cli::cli_abort("Found a NaN in {.arg A}.")
  ## Eliminate explicit zeros (matches Python `A.eliminate_zeros()`).
  A <- Matrix::drop0(A)
  res <- .solve_internal(A, b, c, cone_dict,
                         warm_start = warm_start,
                         solve_method = solve_method,
                         P = P, ...)
  list(x = res$x, y = res$y, s = res$s, info = res$info)
}

#' Solve a cone program and return forward / adjoint derivative operators
#'
#' Mirrors `diffcp.cone_program.solve_and_derivative`. Solves
#'   minimise  c^T x        s.t.  A x + s = b,  s in K
#' (with optional QP `0.5 x^T P x` term) and returns the optimal
#' `(x, y, s)` together with closures `D` (forward) and `DT` (adjoint)
#' that map perturbations of `(A, b, c, [P])` to perturbations of
#' `(x, y, s)` and vice versa.
#'
#' @inheritParams solve_only
#' @param mode Differentiation mode: `"lsqr"` (default), `"dense"`,
#'   or `"lpgd"`.
#' @returns A list with elements `x`, `y`, `s`, `info`, `D`, `DT`.
#'   `info` is the solver-status block returned by the underlying
#'   solver (`status`, `iter`, `solveTime`, `pobj`, ...); see
#'   `solve_only` for the same shape.
#' @export
solve_and_derivative <- function(A, b, c, cone_dict,
                                 P = NULL,
                                 solve_method = "Clarabel",
                                 mode = "lsqr",
                                 warm_start = NULL,
                                 ...) {
  if (!mode %in% c("dense", "lsqr")) {
    cli::cli_abort("Unsupported {.arg mode} = {.val {mode}}; supported: 'dense', 'lsqr'.")
  }
  if (!is.null(P)) {
    ## Mirror Python diffcp/cone_program.py:637-639.  Python diffcp
    ## supports QP only via lpgd modes (lpgd / lpgd_left / lpgd_right);
    ## we have not ported lpgd, so the path errors with the same text.
    cli::cli_abort(c(
      "Dense, lsqr, and lsmr modes currently do not support quadratic objectives.",
      "i" = "Consider switching to {.code mode = 'lpgd'} mode.",
      "x" = "{.code mode = 'lpgd'} is not yet ported in this R version of diffcp."
    ))
  }
  if (any(is.nan(A@x))) cli::cli_abort("Found a NaN in {.arg A}.")
  A <- Matrix::drop0(A)

  res <- .solve_internal(A, b, c, cone_dict,
                         warm_start = warm_start,
                         solve_method = solve_method,
                         P = P, ...)
  x <- res$x; y <- res$y; s <- res$s

  cones <- parse_cone_dict(cone_dict)
  closures <- .make_derivative_closures(A, b, c, x, y, s, cones, mode)

  list(x = x, y = y, s = s, info = res$info,
       D = closures$D, DT = closures$DT)
}
