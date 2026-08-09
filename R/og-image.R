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

OUT <- "."          # where the PNGs land

# palette (styles.scss) --------------------------------------------------------
dark  <- "#141A1C"  # $theme-dark-gray   — page background
lgray <- "#D5D5D5"  # $theme-light-gray  — body text
blue  <- "#60D1F2"  # $theme-light-blue
red   <- "#E76C67"  # $theme-light-red
green <- "#6DCB60"  # $theme-soft-green
gold  <- "#EFC24B"  # $theme-soft-gold

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
plot_variety <- function(pts, poly = NULL, cols = c(gold, green, blue),
                         size = 0.5, alpha = 0.8, curve = FALSE, xlim, ylim) {
  p <- ggplot(pts, aes(x, y, colour = d)) +
    geom_point(size = size, alpha = alpha) +
    scale_colour_gradientn(colours = cols)
  if (curve && !is.null(poly)) {
    p <- p + geom_variety(poly = poly, colour = lgray,
                          linewidth = 0.3, inherit.aes = FALSE)
  }
  p + coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) + bare
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


# render -----------------------------------------------------------------------

set.seed(42)

save_og(
  plot_variety(sample_variety(lemniscate, sd = 0.03, w = 2), poly = lemniscate,
               xlim = c(-1.85, 1.85), ylim = c(-0.971, 0.971)),
  "og-a-variety.png"
)

save_og(plot_streams(field_spiral), "og-b-streams.png")

save_og(
  plot_hdr(rmix(), xlim = c(-3.6, 3.4), ylim = c(-1.7, 1.65)),
  "og-c-hdr.png"
)
