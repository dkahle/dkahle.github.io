# banner.R ---------------------------------------------------------------------
#
# LinkedIn / social banner: sample paths from a Gaussian bridge whose covariance
# slides continuously from Brownian to squared-exponential across its width. Each
# path leaves 0, is jagged early, smooths out, and returns to 0.
#
# Target spec: 1584 x 396 px (4:1), LinkedIn's banner size. Colors come from
# styles.scss so it matches the site.
#
#   Rscript R/banner.R
#
# Needs ggplot2 only.
#
# The construction ------------------------------------------------------------
#
# Two bridge kernels on [0, 1]:
#
#   rough   k_B(s,t) = min(s,t) - s t          Brownian bridge; Holder-1/2, so
#                                              nowhere differentiable — the jagged end.
#   smooth  k_S(s,t) = exp(-(s-t)^2 / 2 l^2)   squared exponential, then conditioned
#                                              to vanish at both endpoints.
#
# NOTE ON THE NAME: the smooth one is the *squared* exponential. "Double
# exponential" usually names the Laplace kernel exp(-|s-t|/l), which is NOT
# differentiable — it is exactly as rough as Brownian motion, so building toward
# it would produce no visible smoothing at all. exp of a *squared* argument is
# what buys infinite differentiability.
#
# Interpolate with a weight w(t) running 1 -> 0:
#
#   k(s,t) = w(s)w(t) k_B(s,t) + (1-w(s))(1-w(t)) k_S(s,t)
#
# This is a genuine non-stationary covariance, not a visual trick: each term is
# D K D for a diagonal D, so each is positive semi-definite, and a sum of PSD
# matrices is PSD. Because both components vanish at t = 0 and t = 1, so does
# the blend — the result is still exactly a bridge. Local roughness at t is
# inherited from whichever component dominates there, which is the effect wanted:
# Brownian on the left, analytic on the right, continuous in between.
#
# Two practical notes:
#
#   * SE kernels are famously ill-conditioned — this one runs ~1e11 — so the
#     Cholesky needs jitter on the diagonal or it fails outright.
#
#   * The two components have very different natural scales (the Brownian bridge
#     peaks at variance 1/4, the SE bridge near 1), so without renormalizing, the
#     right-hand side would tower over the left. Rescaling the covariance to a
#     chosen variance envelope, K -> D K D with D = diag(sqrt(env/diag(K))),
#     fixes the amplitude while leaving the correlation structure — and so the
#     roughness — untouched.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages(library(ggplot2))

OUT  <- "."
SAVE <- !interactive()

dark  <- "#141A1C"  # $theme-dark-gray
green <- "#6DCB60"  # $theme-soft-green

#' Condition a GP to be exactly 0 at the given grid indices.
pin <- function(K, idx, jit = 1e-8) {
  Kuu <- K[idx, idx, drop = FALSE]
  Kfu <- K[, idx, drop = FALSE]
  K - Kfu %*% solve(Kuu + diag(jit, length(idx))) %*% t(Kfu)
}

#' Covariance of the Brownian-to-smooth bridge.
#'
#' @param ell   lengthscale of the smooth component; larger = lazier sweeps.
#' @param shape exponent on w(t) = (1-t)^shape. Larger pushes the handover to the
#'              smooth regime earlier, leaving more of the width elegant and less
#'              of it jagged.
#' @param env_p exponent on the variance envelope (4t(1-t))^env_p. Lower is
#'              flatter, i.e. the paths reach full amplitude sooner.
bridge_kernel <- function(t, ell = 0.20, shape = 1.4, env_p = 0.7) {
  n <- length(t)
  k_rough  <- outer(t, t, pmin) - outer(t, t, "*")
  k_smooth <- pin(exp(-outer(t, t, "-")^2 / (2 * ell^2)), c(1, n))

  w <- (1 - t)^shape
  K <- (w %o% w) * k_rough + ((1 - w) %o% (1 - w)) * k_smooth

  env <- (4 * t * (1 - t))^env_p                  # common amplitude envelope
  s <- sqrt(pmax(env, 0) / pmax(diag(K), 1e-12))
  s[!is.finite(s)] <- 0
  (s %o% s) * K
}

#' Draw m paths. Jitter is not optional here — see the note up top.
sample_paths <- function(K, m, jit = 1e-8) {
  n <- nrow(K)
  L <- tryCatch(chol(K + diag(jit, n)),
                error = function(e) chol(K + diag(jit * 100, n)))
  crossprod(L, matrix(rnorm(n * m), n, m))
}

#' @param glow draw each path three times — wide and faint, then narrower and
#'   brighter — which fakes a bloom without any compositing.
banner <- function(K, t, m = 12, seed = 11, glow = TRUE,
                   col = green, lw = 0.45, alpha = 0.8) {
  set.seed(seed)
  P  <- sample_paths(K, m)
  df <- data.frame(t = rep(t, m), y = as.vector(P),
                   id = factor(rep(seq_len(m), each = length(t))))

  g <- ggplot(df, aes(t, y, group = id))
  if (glow) g <- g +
    geom_line(color = col, linewidth = lw * 6.0, alpha = alpha * 0.06) +
    geom_line(color = col, linewidth = lw * 2.5, alpha = alpha * 0.15)

  g + geom_line(color = col, linewidth = lw, alpha = alpha) +
    coord_cartesian(expand = FALSE) +
    theme_void() +
    theme(plot.background  = element_rect(fill = dark, color = NA),
          panel.background = element_rect(fill = dark, color = NA),
          legend.position  = "none",
          plot.margin      = margin(0, 0, 0, 0))
}

save_banner <- function(p, file, width = 1584, height = 396, dpi = 150) {
  path <- file.path(OUT, file)
  ggsave(path, p, width = width, height = height, units = "px", dpi = dpi, bg = dark)
  cat(sprintf("%-28s %s\n", file, format(file.size(path), big.mark = ",")))
  invisible(path)
}


# 1. data ----------------------------------------------------------------------

n_grid <- 700
tt <- seq(0, 1, length.out = n_grid)
K  <- bridge_kernel(tt, ell = 0.20, shape = 1.4)


# 2. plot ----------------------------------------------------------------------
# LinkedIn overlays the profile photo on the lower-left of the banner — roughly a
# 152px circle centred near x = 8%, hanging off the bottom edge. It clips an
# empty corner at these settings, but re-check with a mock circle if you retune
# shape or env_p, since both move where the dense jagged region sits.

p_banner <- banner(K, tt, m = 12, seed = 11)


# 3. save ----------------------------------------------------------------------

if (SAVE) save_banner(p_banner, "banner.png")
