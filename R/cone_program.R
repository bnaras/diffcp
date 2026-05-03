## DIFFCP SOURCE: diffcp/cone_program.py
##
## Forward solve dispatch (Clarabel / SCS) plus the public
## `solve_only` and `solve_and_derivative` entry points. Derivative
## computation (M operator, LSQR, dense modes) lands in subsequent
## commits; for now `solve_and_derivative` errors with a clear
## phase-pointer message, while `solve_only` is fully functional.

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
    if (any(cone_dict[["s"]] %||% integer(0) > 0L)) {
      cli::cli_abort(c(
        "PSD cones not yet supported via Clarabel in {.pkg diffcp}.",
        "i" = "Clarabel uses upper-triangular ordering; permutation arrives in Phase 2b."
      ))
    }

    ## Clarabel uses the same SCS-style cones list. Translate `f`->`z`
    ## (free cone alias) and drop `ed` (Clarabel does not implement
    ## dual exponential cones natively).
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
    if (!is.null(cone_dict[["ep"]]) && cone_dict[["ep"]] > 0L)
      cones_arg$ep <- as.integer(cone_dict[["ep"]])

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

    if (scs_status != "Solved") {
      if (raise_on_error) {
        cli::cli_abort("Clarabel returned status {.val {status_str}}.")
      } else {
        return(list(x = sol$x, y = sol$z, s = sol$s,
                    info = list(status = scs_status), D = NULL, DT = NULL))
      }
    }

    return(list(
      x = as.numeric(sol$x),
      y = as.numeric(sol$z),
      s = as.numeric(sol$s),
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
    P_csc <- if (is.null(P)) NULL else methods::as(P, "CsparseMatrix")
    sol <- scs::scs(A = A_csc, b = b, obj = c, P = P_csc,
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
  list(x = res$x, y = res$y, s = res$s)
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
#' @returns A list with elements `x`, `y`, `s`, `D`, `DT`.
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
    cli::cli_abort("QP objective {.arg P} not yet supported (Phase 2-QP).")
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

  list(x = x, y = y, s = s, D = closures$D, DT = closures$DT)
}
