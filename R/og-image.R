# og-image.R ------------------------------------------------------------------
#
# Generates candidate Open Graph / social-preview images for kahle.io.
# Target spec: 1200 x 630 px (1.91:1), < ~1 MB, legible when scaled to a
# thumbnail. Colours come from styles.scss so the card matches the site.
#
# Run it from anywhere, including the project root:
#
#   Rscript R/og-image.R
#
# Needs ggplot2, ggvfields, ggdensity, vnorm, and mpoly on the library path.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggvfields)
  library(ggdensity)
  library(vnorm)
  library(mpoly)     # masks ggplot2::vars — harmless, nothing here uses vars()
})

OUT  <- "."              # where the PNGs land
SAVE <- !interactive()   # Rscript writes the PNGs; source()-ing in a session
                         # just builds p_a/p_b/p_c. Override by setting it.

# palette (styles.scss) --------------------------------------------------------
dark  <- "#141A1C"  # $theme-dark-gray   — page background
lgray <- "#D5D5D5"  # $theme-light-gray  — body text
blue  <- "#60D1F2"  # $theme-light-blue
red   <- "#E76C67"  # $theme-light-red
green <- "#6DCB60"  # $theme-soft-green
gold  <- "#EFC24B"  # $theme-soft-gold

# Density ramps for the variety sample, low -> high. Which one you want depends
# on n: warm is fine when the cloud is dense enough to read as a tube, but with
# only ~100 draws it makes the *farthest* points the brightest things on the
# canvas and the sample looks like scatter. fade sinks the tails into the
# background so the eye follows the curve.
dens_warm <- c(gold, green, blue)
dens_fade <- c("#2E3A40", "#4E7F86", blue)

# strip all chrome — the card is pure image, no axes/legend/margins
bare <- theme_void() + theme(
  plot.background  = element_rect(fill = dark, colour = NA),
  panel.background = element_rect(fill = dark, colour = NA),
  legend.position  = "none",
  plot.margin      = margin(0, 0, 0, 0)
)

#' Write a plot at exact OG dimensions.
save_og <- function(p, file, width = 1200, height = 630, dpi = 150) {
  path <- file.path(OUT, file)
  ggsave(path, p, width = width, height = height,
         units = "px", dpi = dpi, bg = dark)
  cat(sprintf("%-28s %s\n", file, format(file.size(path), big.mark = ",")))
  invisible(path)
}





HMC = function (U, grad_U, epsilon, L, current_q) {
  q = current_q
  p = rnorm(length(q), 0, 1)  # independent standard normal variates
  current_p = p
  
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  
  # Alternate full steps for position and momentum
  for (i in 1:L) {
    # Make a full step for the position
    q = q + epsilon * p
    
    # Make a full step for the momentum, except at end of trajectory
    if (i != L) p = p - epsilon * grad_U(q)
  }
  
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q) / 2
  
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (log(runif(1)) < current_U - proposed_U + current_K - proposed_K) {
    return (q)  # accept
  } else {
    return (current_q)  # reject
  }
}


# A. algebraic variety ---------------------------------------------------------
#
# Draws from the variety normal distribution with vnorm: points concentrated on
# the real variety {p = 0} with Gaussian falloff perpendicular to it, coloured
# by the density that actually generated them. So the colour is meaningful —
# blue is the high-density core, gold are genuine tail draws.
#
# Two things that will bite you:
#
#   * rejection = TRUE. The default routes through Stan, which for these curves
#     is far slower than it's worth (it hung a session for 4+ minutes). The
#     rejection sampler does 9k points in ~4 s.
#
#   * w is a BOX WINDOW and it defaults to (-1, 1) on every variable. Leave it
#     and the sample is silently clipped — the lemniscate below runs out to
#     x = ±1.79, so the default window quietly lops off both lobes and you get
#     a plausible-looking but wrong picture. Always set w past the curve's
#     real extent.

#' Sample from the variety normal distribution and score by its density.
#'
#' @param poly an mpoly object, e.g. mp("y^2 - x^3 + x").
#' @param sd   perpendicular spread around the variety. Smaller = tighter tube.
#' @param w    half-width of the sampling box; must exceed the curve's extent.
sample_variety <- function(poly, n = 9000, sd = 0.03, w = 2) {
  s <- rvnorm(n, poly, sd = sd, rejection = TRUE, w = w)
  df <- as.data.frame(s)
  names(df) <- c("x", "y")
  df$d <- pdvnorm(s, poly, sd = sd)
  df
}

#' @param curve overlay the exact variety with vnorm::geom_variety().
#' @param path  optional leapfrog trajectory from hmc_path(), drawn with a
#'              closed arrowhead at the proposal end.
plot_variety <- function(pts, poly = NULL, cols = dens_warm,
                         size = 0.5, alpha = 0.8, curve = FALSE,
                         path = NULL, path_colour = red, path_width = 0.4,
                         head = 0.07, xlim, ylim) {
  p <- ggplot(pts, aes(x, y, colour = d)) +
    geom_point(size = size, alpha = alpha) +
    scale_colour_gradientn(colours = cols)
  if (curve && !is.null(poly)) {
    p <- p + geom_variety(poly = poly, colour = lgray,
                          linewidth = 0.3, inherit.aes = FALSE)
  }
  if (!is.null(path)) {
    p <- p + geom_path(
      data = path, aes(x, y), inherit.aes = FALSE,
      colour = path_colour, linewidth = path_width,
      arrow = arrow(type = "closed", angle = 22, length = unit(head, "in"))
    )
  }
  p + coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) + bare
}


# One HMC transition, kept as a path -------------------------------------------
#
# HMC() above returns only the endpoint, so there's nothing to draw. hmc_path()
# is the same leapfrog integrator recording q at every step, which is the thing
# that visibly weaves across the variety.
#
# The potential is the variety normal's: U = g^2 / (2 sd^2), so grad U =
# g * grad(g) / sd^2. Derivatives come from mpoly symbolically — exact, which
# matters because leapfrog on a stiff potential is unforgiving.
#
# Tuning, the short version:
#
#   * WEAVE AMPLITUDE IS ROUGHLY sd/|grad g|, so it is set by sd, not by the
#     integrator. At sd = 0.03 the oscillation is ~0.006 units — invisible. The
#     weave only reads once the tube is wide enough, hence sd = 0.12 below. The
#     sample and the path must share sd or the path wanders outside the cloud.
#
#   * epsilon has a hard stability ceiling near 2 sd/|grad g|. Above it the
#     trajectory blows up to NaN and the acceptance test silently returns NA
#     (that is what isTRUE() below is guarding). At sd = 0.03, eps = 0.010
#     already diverges.

#' Build the variety normal potential and its gradient for an mpoly.
variety_potential <- function(poly, sd, vars = c("x", "y")) {
  g  <- as.function(poly, varorder = vars, silent = TRUE)
  gx <- as.function(deriv(poly, vars[1]), varorder = vars, silent = TRUE)
  gy <- as.function(deriv(poly, vars[2]), varorder = vars, silent = TRUE)
  list(
    U      = function(q) g(q)^2 / (2 * sd^2),
    grad_U = function(q) g(q) * c(gx(q), gy(q)) / sd^2
  )
}

#' One HMC transition, returning the whole leapfrog trajectory.
#'
#' @return data frame of the path, with attributes "accept" and "arclen".
hmc_path <- function(poly, sd, epsilon, L, start, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  pot <- variety_potential(poly, sd)
  U <- pot$U; grad_U <- pot$grad_U

  q <- start
  p <- rnorm(length(q), 0, 1)
  current_q <- q; current_p <- p

  path <- matrix(NA_real_, L + 1, 2)
  path[1, ] <- q

  p <- p - epsilon * grad_U(q) / 2
  for (i in 1:L) {
    q <- q + epsilon * p
    path[i + 1, ] <- q
    if (i != L) p <- p - epsilon * grad_U(q)
  }
  p <- -(p - epsilon * grad_U(q) / 2)

  accept <- isTRUE(
    log(runif(1)) < U(current_q) - U(q) + sum(current_p^2) / 2 - sum(p^2) / 2
  )

  df <- as.data.frame(path)
  names(df) <- c("x", "y")
  if (anyNA(df)) warning("leapfrog diverged — lower epsilon")
  attr(df, "accept") <- accept
  attr(df, "arclen") <- sum(sqrt(diff(df$x)^2 + diff(df$y)^2))
  df
}

#' A point exactly on the lemniscate r^2 = a cos(2 theta), to start the chain.
on_lemniscate <- function(theta, a = 3.2) {
  r <- sqrt(a * cos(2 * theta))
  c(r * cos(theta), r * sin(theta))
}

# Curves to try. The frame is 1.91:1, so wide curves fit best; set w and the
# coord_fixed window to match the curve's real extent (printed by sample_variety
# if you check range()).
lemniscate <- mp("x^4 + 2 x^2 y^2 + y^4 - 3.2 x^2 + 3.2 y^2")  # x to ±1.79
elliptic   <- mp("y^2 - x^3 + x")                # oval + branch; lots of dead space
trifolium  <- mp("x^4 + 2 x^2 y^2 + y^4 - x^3 + 3 x y^2")


# B. ggvfields stream field ----------------------------------------------------
#
# T is the integration time and the whole ballgame: the default leaves stubs
# that read as an arrow grid, T = 12 packs the frame into a solid vortex.
# T ≈ 1.5 is the readable middle. arrow = NULL drops the arrowheads, which
# otherwise dominate at this scale.

plot_streams <- function(fun, xlim = c(-4.4, 4.4), ylim = c(-2.4, 2.4),
                         n = 15, T = 1.5, colour = blue,
                         linewidth = 0.55, alpha = 0.8,
                         view_x = c(-3.8, 3.8), view_y = c(-2, 2)) {
  ggplot() +
    geom_stream_field(
      fun = fun, xlim = xlim, ylim = ylim, n = n, T = T,
      normalize = FALSE, arrow = NULL,
      linewidth = linewidth, colour = colour, alpha = alpha
    ) +
    coord_fixed(xlim = view_x, ylim = view_y, expand = FALSE) +
    bare
}

# sample fields — return c(dx, dy) for a point v = c(x, y)
field_spiral <- function(v) c(-v[2] + 0.35 * sin(1.1 * v[1]),
                              0.75 * v[1] + 0.35 * cos(1.1 * v[2]))
field_saddle <- function(v) c(v[1] + 0.4 * v[2], -v[2] + 0.4 * sin(v[1]))
field_dipole <- function(v) c(v[1]^2 - v[2]^2 - 1, 2 * v[1] * v[2])


# C. ggdensity HDRs ------------------------------------------------------------

plot_hdr <- function(df, probs = c(.99, .95, .8, .5), fill = blue,
                     points = TRUE, xlim, ylim) {
  p <- ggplot(df, aes(x, y)) + geom_hdr(probs = probs, fill = fill)
  if (points) p <- p + geom_point(size = 0.18, alpha = 0.25, colour = lgray)
  p + coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) + bare
}

#' Draw from a bivariate normal mixture.
rmix <- function(n = 3000,
                 mu_x = c(-1.5, 1.4, 0.2), mu_y = c(-0.35, 0.45, -0.6),
                 sd_x = c(.75, .6, .45),   sd_y = c(.45, .4, .3),
                 w    = c(.40, .35, .25)) {
  k <- sample(length(w), n, replace = TRUE, prob = w)
  data.frame(x = rnorm(n, mu_x[k], sd_x[k]),
             y = rnorm(n, mu_y[k], sd_y[k]))
}


# 1. data ----------------------------------------------------------------------
# The slow step: rvnorm() runs a few seconds. Do this once, then iterate on the
# plots below without paying for it again.

set.seed(42)

# sd is shared by the sample and the HMC path so the trajectory weaves through
# the same band the points occupy — see the note above hmc_path().
SD_A <- 0.12

d_variety <- sample_variety(lemniscate, n = 100, sd = SD_A, w = 2)

# seed 7 is chosen, not lucky: it traverses the whole lemniscate — left lobe,
# through the node, around the right — instead of stalling in one lobe.
d_path <- hmc_path(lemniscate, sd = SD_A, epsilon = 0.010, L = 340,
                   start = on_lemniscate(0.45), seed = 7)

d_hdr <- rmix()

# The stream field has no data step — geom_stream_field() integrates the field
# at draw time, so its "data" is the function field_spiral itself.


# 2. plots ---------------------------------------------------------------------
# Cheap to rebuild. Tweak here and print p_a / p_b / p_c to look at one.

p_a <- plot_variety(d_variety, poly = lemniscate, curve = TRUE,
                    path = d_path, size = 0.8, cols = dens_fade,
                    xlim = c(-1.85, 1.85), ylim = c(-0.971, 0.971))

p_b <- plot_streams(field_spiral)

p_c <- plot_hdr(d_hdr, xlim = c(-3.6, 3.4), ylim = c(-1.7, 1.65))


# 3. save ----------------------------------------------------------------------
# A plot printed to the IDE pane is NOT the card. These all use coord_fixed(),
# so the pane letterboxes them at whatever aspect the pane happens to be —
# composition only reads true at 1200x630. Check the PNG before settling on a
# design.

if (SAVE) {
  save_og(p_a, "og-a-variety.png")
  save_og(p_b, "og-b-streams.png")
  save_og(p_c, "og-c-hdr.png")
}
