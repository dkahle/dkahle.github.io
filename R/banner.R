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
#
# HOW THE SETTINGS WERE ARRIVED AT ---------------------------------------------
#
# Final: ell 0.20, share (1-t)^4, amp_exp(0.25), m = 16 paths over 8 colors,
# seed 56, glow off, refined k = 4 (700 -> 2797 points). Recorded because in
# every case the obvious choice was the wrong one:
#
#   share, the transition. Blending the raw kernels under a single weight gave
#     no real control at all (see WHY NORMALISE FIRST above). Within the family,
#     BIGGER exponents smooth sooner — so (1-t) is its slowest member, not its
#     fastest, which is the reverse of what "linear" suggests. Ramping share to
#     exactly zero partway across scores best of anything tried (roughness 0.004
#     past the cutoff) and looks worst: beyond the cutoff the process is purely
#     squared-exponential, whose correlation length is far longer, so paths snap
#     from jagged into clean arcs. (1-t)^4 clears the right ~40% with no seam.
#
#   amp, the envelope. Originally an exponent on the VARIANCE — a trap, since
#     the visible spread is its square root, so a genuinely quadratic variance
#     drew a linear cone. It now takes the amplitude directly. Polynomials of
#     any degree still leave a corner where the paths lift off the baseline;
#     exp(-c/(t(1-t))) is C-infinity, every derivative vanishing at both ends,
#     so there is no corner at any magnification. c = 0.25 holds the envelope
#     near 1% of full height at t = 0.05.
#
#   seed. 90 seeds scored on four things: whether the bulk of paths fills the
#     frame (rather than one spike squashing the rest), spread through the
#     smooth half, spread through the jagged half, and balance about the
#     centerline. 37 passed, which only rules out the obviously bad, so the top
#     12 were rendered and judged by eye. 56 has the fullest body and the
#     cleanest sweeps once smooth.
#
#   grid. 700 points while selecting, then refine() to 2797 by conditional
#     simulation — which preserves the chosen realization instead of redrawing
#     it. See refine() for why resampling would not.
#
# Reproduce exactly with: Rscript R/banner.R
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
#'   the midpoint. BIGGER EXPONENTS SMOOTH SOONER, which is the opposite of what
#'   "linear" suggests — (1-t) is the slowest of the family, not the fastest,
#'   and is still faintly jagged at the right edge because Brownian dominates at
#'   small scales even at a low share. (1-t)^4 clears the right ~40% into smooth
#'   sweeps; ^2 stays rough past the middle; beyond ^6 the handover compresses.
#'
#'   Do NOT ramp share linearly to exactly zero part-way across, e.g.
#'   pmax(0, 1 - t/0.6). It looks like the fastest option on paper and it is
#'   measurably clean — roughness falls to 0.004 past the cutoff — but the
#'   process becomes purely squared-exponential there, whose correlation length
#'   is far longer, so the paths visibly snap from jagged into clean arcs. It
#'   reads as a switch being flipped, not a transition.
#' @param amp function of t giving the AMPLITUDE envelope — the spread of the
#'   paths you actually see, not the variance. It is parameterised this way on
#'   purpose: an earlier version took an exponent on the variance, and since the
#'   visible spread is its square root, "quadratic variance" silently drew a
#'   linear cone. Squaring happens internally now, so what you write is what you
#'   see. Use amp_poly() or amp_exp() below.
bridge_kernel <- function(t, ell = 0.20, share = function(t) (1 - t)^4,
                          amp = amp_exp(0.25)) {
  n <- length(t)
  C_rough  <- to_corr(outer(t, t, pmin) - outer(t, t, "*"))
  C_smooth <- to_corr(pin(exp(-outer(t, t, "-")^2 / (2 * ell^2)), c(1, n)))

  a <- sqrt(pmax(pmin(share(t), 1), 0))
  b <- sqrt(1 - a^2)
  C <- (a %o% a) * C_rough + (b %o% b) * C_smooth

  e <- amp(t)                                     # amplitude, squared by the outer product
  (e %o% e) * C
}

#' Polynomial amplitude envelope: the paths open as t^p at each end.
amp_poly <- function(p = 2) function(t) (4 * t * (1 - t))^p

#' Exponential amplitude envelope — exp(-c / (t(1-t))), normalised to peak at 1.
#'
#' This is the C-infinity bump: every derivative vanishes at both endpoints, so
#' it is flatter near the ends than any polynomial, however high the degree. The
#' paths lift off the baseline with no visible corner at all, then swell into a
#' fuller body than a polynomial gives, because once clear of the ends the bump
#' rises quickly and plateaus.
#'
#' c controls how long it hugs zero. At c = 0.10 the envelope is already at 18%
#' of full height by t = 0.05; at 0.25 it is 1%; at 0.35, 0.3%. Larger means a
#' longer flat run and a more abrupt swell.
amp_exp <- function(c = 0.25) function(t) {
  z <- exp(-c / (t * (1 - t)))
  z[!is.finite(z)] <- 0
  z / max(z)
}

#' Refine an existing realization onto a finer grid, keeping it unchanged.
#'
#' Resampling on a finer grid does NOT give you the same picture at higher
#' resolution — sample_paths() draws n x m normals, so changing n changes the
#' draw and you get a different realization from the same seed. To add detail to
#' a realization you already like, condition on it: hold X fixed at the coarse
#' nodes and simulate the new points from the conditional law.
#'
#' Worth doing because the rough half is Brownian, which has structure at every
#' scale — a finer grid is not cosmetic there, it resolves texture that a coarse
#' grid simply never drew. The smooth half is unaffected, as it should be.
#'
#' @param k points inserted between each pair of coarse nodes; k = 4 takes 700
#'   nodes to 2797.
#' @return list(t, X, idx); X[idx, ] reproduces X_c to within the jitter.
refine <- function(t_c, X_c, k = 4, ell = 0.20,
                   share = function(t) (1 - t)^4, amp = amp_exp(0.25),
                   jit = 1e-9) {
  t_f <- sort(unique(c(t_c, as.vector(sapply(seq_len(length(t_c) - 1), function(i)
           seq(t_c[i], t_c[i + 1], length.out = k + 1)[-c(1, k + 1)])))))
  idx <- match(t_c, t_f)
  stopifnot(!anyNA(idx))

  K   <- bridge_kernel(t_f, ell = ell, share = share, amp = amp)
  Kcc <- K[idx, idx]; Kfc <- K[, idx]
  A   <- Kfc %*% solve(Kcc + diag(jit, length(idx)))
  mu  <- A %*% X_c
  Sig <- K - A %*% t(Kfc)

  # Sig is singular by construction — conditional variance is exactly 0 at every
  # coarse node — so a plain Cholesky fails. Pivoting handles the rank deficit.
  R <- suppressWarnings(chol(Sig + diag(jit, nrow(Sig)), pivot = TRUE))
  p <- attr(R, "pivot")
  y <- matrix(0, nrow(Sig), ncol(X_c))
  y[p, ] <- crossprod(R, matrix(rnorm(nrow(Sig) * ncol(X_c)), nrow(Sig)))

  list(t = t_f, X = mu + y, idx = idx)
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
#' @param paths optional precomputed path matrix (e.g. from refine()); when
#'   given, K and seed are ignored and t must match its rows.
banner <- function(K, t, m = 16, seed = 16, glow = TRUE, cols = soft8,
                   lw = 0.45, alpha = 0.85, paths = NULL) {
  if (is.null(paths)) { set.seed(seed); paths <- sample_paths(K, m) }
  P <- paths
  m <- ncol(P)
  stopifnot(nrow(P) == length(t))
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
K  <- bridge_kernel(tt, ell = 0.20, share = function(t) (1 - t)^4,
                    amp = amp_exp(0.25))


# 2. plot ----------------------------------------------------------------------
# LinkedIn overlays the profile photo on the lower-left of the banner — roughly a
# 152px circle centred near x = 8%, hanging off the bottom edge. It clips an
# empty corner at these settings, but re-check with a mock circle if you retune
# share or amp, since both move where the dense jagged region sits.

SEED <- 56
set.seed(SEED); X_coarse <- sample_paths(K, 16)
set.seed(SEED); fine <- refine(tt, X_coarse, k = 4)   # 700 -> 2797 points

(p_banner <- banner(K = NULL, fine$t, paths = fine$X, glow = FALSE))


# 3. save ----------------------------------------------------------------------

if (SAVE) save_banner(p_banner, "banner.png")
