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
# Normalise each to unit variance first, then interpolate the CORRELATIONS with
# a(t)^2 + b(t)^2 = 1:
#
#   k(s,t) = a(s)a(t) C_B(s,t) + b(s)b(t) C_S(s,t),   then scale to an envelope
#
# This is a genuine non-stationary covariance, not a visual trick: each term is
# D C D for a diagonal D, so each is positive semi-definite, and a sum of PSD
# matrices is PSD. Because the envelope vanishes at t = 0 and t = 1, so does the
# blend — the result is still exactly a bridge. Local roughness at t is inherited
# from whichever component dominates there: Brownian on the left, analytic on the
# right, continuous in between.
#
# WHY NORMALISE FIRST. Blending the raw kernels instead, with a single weight
# w(t), does not let you steer the transition. The rough component's contribution
# is w(t)^2 k_B(t,t), and k_B(t,t) = t(1-t) is itself already collapsing toward
# t = 1, so the two decays compound and the handover finishes far earlier than
# w(t) suggests — measurably, roughness was down to a quarter of its peak by
# t = 0.4 and essentially gone by 0.6. Normalising to correlations first strips
# that hidden factor out, so a(t)^2 IS the share of variance that is Brownian at
# t, and the transition goes exactly where you put it.
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

# The four site accents from styles.scss, plus four more picked to sit at the
# same lightness and saturation so no single path jumps forward on the dark
# ground. Paths cycle through these.
soft8 <- c("#60D1F2",  # $theme-light-blue
           "#6DCB60",  # $theme-soft-green
           "#EFC24B",  # $theme-soft-gold
           "#E76C67",  # $theme-light-red
           "#5BD6B0",  # aqua
           "#9B8CE8",  # violet
           "#F0A35E",  # orange
           "#E884B8")  # pink

#' Condition a GP to be exactly 0 at the given grid indices.
pin <- function(K, idx, jit = 1e-8) {
  Kuu <- K[idx, idx, drop = FALSE]
  Kfu <- K[, idx, drop = FALSE]
  K - Kfu %*% solve(Kuu + diag(jit, length(idx))) %*% t(Kfu)
}

#' Rescale a kernel to unit variance, i.e. take its correlation matrix.
to_corr <- function(K) {
  d <- sqrt(diag(K)); d[d < 1e-9] <- 1e-9
  C <- K / outer(d, d)
  C[!is.finite(C)] <- 0
  diag(C) <- 1
  C
}

#' Covariance of the Brownian-to-smooth bridge.
#'
#' @param ell   lengthscale of the smooth component; larger = lazier sweeps.
#' @param share function of t giving the fraction of variance that is Brownian.
#'   Must run 1 -> 0. This is the transition control, and it means what it says:
#'   share(0.5) = 0.25 puts a quarter of the variance in the rough component at
#'   the midpoint. Bigger exponents finish the handover sooner. (1-t)^2 keeps
#'   visible roughness past the middle and still arrives smooth; (1-t) is so
#'   gradual it is still faintly jagged at the right edge, because Brownian
#'   dominates at small scales even at a low share.
#' @param env_p exponent on the variance envelope (4t(1-t))^env_p. Lower is
#'              flatter, i.e. the paths reach full amplitude sooner.
bridge_kernel <- function(t, ell = 0.20, share = function(t) (1 - t)^2,
                          env_p = 0.7) {
  n <- length(t)
  C_rough  <- to_corr(outer(t, t, pmin) - outer(t, t, "*"))
  C_smooth <- to_corr(pin(exp(-outer(t, t, "-")^2 / (2 * ell^2)), c(1, n)))

  a <- sqrt(pmax(pmin(share(t), 1), 0))
  b <- sqrt(1 - a^2)
  C <- (a %o% a) * C_rough + (b %o% b) * C_smooth

  e <- sqrt((4 * t * (1 - t))^env_p)              # common amplitude envelope
  (e %o% e) * C
}

#' Draw m paths. Jitter is not optional here — see the note up top.
sample_paths <- function(K, m, jit = 1e-8) {
  n <- nrow(K)
  L <- tryCatch(chol(K + diag(jit, n)),
                error = function(e) chol(K + diag(jit * 100, n)))
  crossprod(L, matrix(rnorm(n * m), n, m))
}

#' @param cols colors to cycle the paths through; with m > length(cols) each
#'   color carries several paths.
#' @param glow draw each path three times — wide and faint, then narrower and
#'   brighter — which fakes a bloom without any compositing. It earns its keep
#'   more with several colors than with one, since the halo separates paths
#'   where they cross.
banner <- function(K, t, m = 16, seed = 16, glow = TRUE, cols = soft8,
                   lw = 0.45, alpha = 0.85) {
  set.seed(seed)
  P  <- sample_paths(K, m)
  df <- data.frame(t = rep(t, m), y = as.vector(P),
                   id  = factor(rep(seq_len(m), each = length(t))),
                   col = rep(rep(cols, length.out = m), each = length(t)))

  g <- ggplot(df, aes(t, y, group = id, color = col))
  if (glow) g <- g +
    geom_line(linewidth = lw * 6.0, alpha = alpha * 0.05) +
    geom_line(linewidth = lw * 2.5, alpha = alpha * 0.13)

  g + geom_line(linewidth = lw, alpha = alpha) +
    scale_color_identity() +
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
K  <- bridge_kernel(tt, ell = 0.20, share = function(t) (1 - t)^2)


# 2. plot ----------------------------------------------------------------------
# LinkedIn overlays the profile photo on the lower-left of the banner — roughly a
# 152px circle centred near x = 8%, hanging off the bottom edge. It clips an
# empty corner at these settings, but re-check with a mock circle if you retune
# shape or env_p, since both move where the dense jagged region sits.

(p_banner <- banner(K, tt, m = 16, seed = 16))


# 3. save ----------------------------------------------------------------------

if (SAVE) save_banner(p_banner, "banner.png")
