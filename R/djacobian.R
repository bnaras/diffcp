## DIFFCP SOURCE: diffcp/cpp/src/dprojection.cpp @ v1.1.8 (6a78143) (C++ extension).
## Upstream pin: see inst/UPSTREAM.dcf / upstream_info().
##
## Per-cone Jacobians of the projection operator.
##
## For each cone K, dpi_K(x) is the (sub-)derivative of Pi_K at x: a
## linear operator from R^n -> R^n. We expose each as a *closure*
## with `matvec(v)` and `rmatvec(v)` methods (the projection's
## Jacobian is symmetric for self-dual cones, so matvec == rmatvec
## there; we still expose both for uniformity with the more general
## non-self-dual cases.)

# -- Zero cone ----------------------------------------------------
## Pi_0(x) = 0 (primal). dPi = 0.
## Pi_{R^n}(x) = x (dual). dPi = I.
.dproj_zero <- function(x, dual = FALSE) {
  n <- length(x)
  if (dual) {
    ## identity
    list(
      matvec  = function(v) v,
      rmatvec = function(v) v,
      n = n
    )
  } else {
    ## zero
    list(
      matvec  = function(v) numeric(n),
      rmatvec = function(v) numeric(n),
      n = n
    )
  }
}

# -- Nonneg orthant -----------------------------------------------
## Pi_{R_+}(x)_i = max(x_i, 0). dPi = diag(1[x > 0]).
## Self-dual: same Jacobian for primal and dual.
.dproj_nonneg <- function(x, dual = FALSE) {
  d <- as.numeric(x > 0)
  list(
    matvec  = function(v) d * v,
    rmatvec = function(v) d * v,
    n = length(x)
  )
}

# -- Single SOC -----------------------------------------------------
## Q^{n+1} = { (t, z) : ||z|| <= t }. Self-dual.
##
## DIFFCP SOURCE: diffcp/cpp/src/dprojection.cpp (SOC branch).
## At x = (t, z), let n = ||z||.
##   - If x is in interior of K (n < t): dPi = I.
##   - If x is in interior of -K (n < -t): dPi = 0.
##   - On the boundary or in the "transition" region, dPi has the
##     closed form below (a low-rank update of a smoothed identity).
##
## The non-trivial branch (norm_z > |t|): the projection writes
##   Pi(x) = (alpha, alpha * z / n)  where alpha = (t + n) / 2.
## Differentiating gives a 2x2 block structure between t and the
## radial direction z/n, plus an isotropic component on the
## tangential complement.
.dproj_soc_one <- function(x, dual = FALSE) {
  n <- length(x)
  if (n == 0L) {
    return(list(matvec = function(v) v, rmatvec = function(v) v, n = 0L))
  }
  t <- x[1L]
  z <- x[-1L]
  norm_z <- sqrt(sum(z * z))

  if (norm_z <= t) {
    ## Interior of K -> identity
    return(list(
      matvec  = function(v) v,
      rmatvec = function(v) v,
      n = n
    ))
  }
  if (norm_z <= -t) {
    ## Interior of -K (= polar of K) -> zero
    return(list(
      matvec  = function(v) numeric(n),
      rmatvec = function(v) numeric(n),
      n = n
    ))
  }

  ## Boundary / transition region.
  ## D = 0.5 * (I + [t/||z|| 0; 0 (t + ||z||)/(2||z||) * I - t/(2||z||) * z z^T / ||z||^2])
  ## Hand-derived from differentiating the closed-form projection;
  ## standard reference is the SCS / diffcp source.
  s <- t / norm_z
  zhat <- z / norm_z
  ## Block matvec:
  ##   v0 := v[1], w := v[-1]
  ##   out0 = 0.5 * (v0 + s * (zhat . w))
  ##   out_w = 0.5 * (s * v0 * zhat + (1 + s) * w - s * (zhat . w) * zhat)
  matvec <- function(v) {
    v0 <- v[1L]
    w  <- v[-1L]
    zw <- sum(zhat * w)
    out0 <- 0.5 * (v0 + zw)
    out_w <- 0.5 * (v0 * zhat + (1 + s) * w - (1 - s) * zw * zhat)
    c(out0, out_w)
  }
  ## SOC projection's Jacobian is symmetric (self-dual cone), so
  ## rmatvec == matvec.
  list(matvec = matvec, rmatvec = matvec, n = n)
}

# -- Cartesian product Jacobian -----------------------------------
## Builds a block-diagonal Jacobian operator from per-cone-block
## Jacobians, returned as a closure with matvec / rmatvec.
##
## Mirrors `diffcp._diffcp.dprojection(v, cones_parsed, dual)` (Python
## C++ extension). Only Zero, Nonneg, and SOC are supported in this
## development version; PSD and EXP arrive in subsequent commits.
.dprojection <- function(x, cones, dual = FALSE) {
  blocks <- list()
  offset <- 0L
  for (entry in cones) {
    cone <- entry$name
    sz   <- entry$sizes
    if (length(sz) == 0L || sum(sz) == 0L) next
    if (cone == "q") {
      for (d in sz) {
        idx <- (offset + 1L):(offset + d)
        blocks[[length(blocks) + 1L]] <- list(
          idx  = idx,
          jac  = .dproj_soc_one(x[idx], dual = dual)
        )
        offset <- offset + d
      }
    } else if (cone == "s" || cone == "ep" || cone == "ed") {
      cli::cli_abort("Jacobian for cone {.val {cone}} not yet implemented.")
    } else {
      for (d in sz) {
        idx <- (offset + 1L):(offset + d)
        jac <- switch(cone,
          z = .dproj_zero(x[idx], dual = dual),
          l = .dproj_nonneg(x[idx], dual = dual),
          cli::cli_abort("Unknown cone {.val {cone}}.")
        )
        blocks[[length(blocks) + 1L]] <- list(idx = idx, jac = jac)
        offset <- offset + d
      }
    }
  }
  total <- offset

  matvec <- function(v) {
    out <- numeric(total)
    for (b in blocks) {
      out[b$idx] <- b$jac$matvec(v[b$idx])
    }
    out
  }
  rmatvec <- function(v) {
    out <- numeric(total)
    for (b in blocks) {
      out[b$idx] <- b$jac$rmatvec(v[b$idx])
    }
    out
  }
  list(matvec = matvec, rmatvec = rmatvec, n = total)
}
