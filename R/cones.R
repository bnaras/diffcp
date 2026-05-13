## DIFFCP SOURCE: diffcp/cones.py @ v1.1.8 (6a78143) (Python source of truth).
## Upstream pin: see inst/UPSTREAM.dcf / upstream_info().
##
## Cone projections: each `proj_*` function takes a numeric vector `x`
## of the appropriate length and returns Pi_K(x), the projection onto
## the cone K. Dual projections Pi_{K^*} are computed via the Moreau
## decomposition Pi_{K^*}(x) = x + Pi_K(-x) for self-dual / non-self-dual
## cases as needed.

#' Parse an SCS-style cone dictionary into an ordered list
#'
#' Mirrors `diffcp.cones.parse_cone_dict`.
#'
#' @param cone_dict A named list with keys among `"z"`, `"l"`, `"q"`,
#'   `"s"`, `"ep"`, `"ed"`. Values are scalar dimensions (for `z`, `l`,
#'   `ep`, `ed`) or integer vectors of dimensions (for `q`, `s`).
#' @returns A list of `list(name, size)` pairs, in the canonical SCS
#'   cone order (`CONES`).
#' @export
parse_cone_dict <- function(cone_dict) {
  out <- list()
  for (cone in CONES) {
    if (!is.null(cone_dict[[cone]])) {
      ## C++ wrapper expects `sizes` to always be an integer vector,
      ## even for scalar-dim cones (z, l, ep, ed). Match that schema
      ## here so R and C++ paths consume the same structure.
      sz <- as.integer(cone_dict[[cone]])
      out[[length(out) + 1L]] <- list(name = cone, sizes = sz)
    }
  }
  out
}

# -- Projection onto the zero cone --------------------------------
## Pi_{0}(x) = 0 (primal); Pi_{R^n}(x) = x (dual: trivial since K^*=R^n).
.proj_zero <- function(x, dual = FALSE) {
  if (dual) x else numeric(length(x))
}

# -- Projection onto the nonneg orthant ---------------------------
## Self-dual: Pi_{R_+}(x) = max(x, 0).
.proj_nonneg <- function(x, dual = FALSE) {
  pmax(x, 0)
}

# -- Projection onto a single second-order (Lorentz) cone ---------
## Q^{n+1} = { (t, z) in R x R^n : ||z||_2 <= t }.
## Self-dual: same projection for primal and dual.
##
## DIFFCP SOURCE: diffcp/cones.py:88-97 (_proj branch for SOC).
.proj_soc_one <- function(x, dual = FALSE) {
  if (length(x) == 0L) return(x)
  t <- x[1L]
  z <- x[-1L]
  norm_z <- sqrt(sum(z * z))
  if (norm_z <= t || isTRUE(all.equal(norm_z, t, tolerance = 1e-8))) {
    return(x)
  }
  if (norm_z <= -t) {
    return(numeric(length(x)))
  }
  scale <- 0.5 * (1 + t / norm_z)
  c(scale * norm_z, scale * z)
}

# -- Projection onto a Cartesian product of SOC cones -------------
## `dims` is an integer vector; the input `x` is split into chunks of
## sizes `dims` and each is projected independently.
.proj_soc <- function(x, dims, dual = FALSE) {
  out <- numeric(length(x))
  offset <- 0L
  for (d in dims) {
    block <- x[(offset + 1L):(offset + d)]
    out[(offset + 1L):(offset + d)] <- .proj_soc_one(block, dual = dual)
    offset <- offset + d
  }
  out
}

# -- Projection onto the PSD cone ---------------------------------
## DIFFCP SOURCE: diffcp/cones.py:99-102 (PSD branch of `_proj`).
##
## PSD is self-dual, so the `dual` argument is ignored.  `x` is the
## SCS-style packed lower triangle (length n*(n+1)/2 with sqrt(2)
## scaling on off-diagonals).  We unpack to a symmetric matrix,
## clip eigenvalues to [0, +Inf), and re-pack.
.proj_psd <- function(x, dim) {
  X <- unvec_symm(x, dim)
  eig <- eigen(X, symmetric = TRUE)
  X_proj <- eig$vectors %*% (pmax(eig$values, 0) * t(eig$vectors))
  vec_symm(X_proj)
}

# -- Single-cone dispatch -----------------------------------------
## Internal helper used by `pi()` to project one block.
.proj_one_block <- function(x, cone, dual = FALSE) {
  switch(cone,
    z  = .proj_zero(x, dual),
    l  = .proj_nonneg(x, dual),
    q  = .proj_soc_one(x, dual),
    ## PSD, EXP, EXP_DUAL deferred to a subsequent commit.
    cli::cli_abort(
      "Cone {.val {cone}} not yet implemented in this development version."
    )
  )
}

#' Projection onto a Cartesian product of cones
#'
#' Mirrors `diffcp.cones.pi`.
#'
#' @param x A numeric vector.
#' @param cones A list of `list(name, size)` pairs, as produced by
#'   `parse_cone_dict`.
#' @param dual If `TRUE`, project onto the dual cones.
#' @returns A numeric vector of the same length as `x`.
#' @export
pi <- function(x, cones, dual = FALSE) {
  out <- numeric(length(x))
  offset <- 0L
  for (entry in cones) {
    cone <- entry$name
    sz   <- entry$sizes
    if (length(sz) == 0L || sum(sz) == 0L) next
    if (cone == "q") {
      ## SOC: each entry is one cone of size `d`.
      for (d in sz) {
        out[(offset + 1L):(offset + d)] <-
          .proj_soc_one(x[(offset + 1L):(offset + d)], dual = dual)
        offset <- offset + d
      }
    } else if (cone == "s") {
      ## PSD: each entry is the matrix dimension; vector size is dim*(dim+1)/2.
      for (d in sz) {
        ## NOTE %/% has higher precedence than * in R; ALWAYS parenthesize.
        sz_vec <- (d * (d + 1L)) %/% 2L
        block_idx <- (offset + 1L):(offset + sz_vec)
        out[block_idx] <- .proj_psd(x[block_idx], d)
        offset <- offset + sz_vec
      }
    } else if (cone == "ep" || cone == "ed") {
      ## DIFFCP SOURCE: diffcp/cones.py:80-114 (EXP / EXP_DUAL branch
      ## of _proj).  EXP_DUAL is handled by toggling `dual`:
      ##   EXP_DUAL, dual=F  ->  EXP, dual=T
      ##   EXP_DUAL, dual=T  ->  EXP, dual=F
      ## Within EXP itself, dual=T uses Moreau:
      ##   Pi_{K_exp^*}(x) = x + Pi_{K_exp}(-x).
      eff_dual <- if (cone == "ed") !dual else dual
      for (count in sz) {
        for (k in seq_len(count)) {
          idx <- (offset + 1L):(offset + 3L)
          xi <- x[idx]
          if (eff_dual) xi <- -xi
          proj <- as.numeric(cpp_project_exp_cone(xi))
          out[idx] <- if (eff_dual) x[idx] + proj else proj
          offset <- offset + 3L
        }
      }
    } else {
      ## Scalar-dim cones: zero, nonneg.
      for (d in sz) {
        block_idx <- (offset + 1L):(offset + d)
        out[block_idx] <- .proj_one_block(x[block_idx], cone, dual = dual)
        offset <- offset + d
      }
    }
  }
  out
}
