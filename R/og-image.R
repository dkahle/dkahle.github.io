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
# Needs ggplot2, ggvfields, and ggdensity on the library path.
# ------------------------------------------------------------------------------

library(ggplot2)
library(ggvfields)
library(ggdensity)

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
# Rejection-sample points near the real variety {g = 0}, then colour by
# distance to it.
#
# The one non-obvious bit: DON'T filter on |g| < tol. The width of that band
# scales like tol/||grad g||, so it balloons wherever the gradient is small —
# which is what made the first attempt blob at the centre. Dividing by the
# gradient norm gives a first-order distance estimate and therefore a tube of
# roughly uniform width. Singular points (where grad g = 0, e.g. the node of a
# lemniscate) will still thicken; that's real, not an artefact.

#' Sample points near the real variety of g.
#'
#' @param g   function(x, y) — the defining polynomial, vectorised.
#' @param tol tube half-width, in units of approximate distance to the curve.
#' @param n   how many points to keep.
#' @param N   how many candidates to draw before rejecting (raise if you get
#'            fewer than n points back).
sample_variety <- function(g, xlim, ylim, tol = 0.045, n = 7000, N = 6e6) {
  x <- runif(N, xlim[1], xlim[2])
  y <- runif(N, ylim[1], ylim[2])

  h  <- 1e-5                                    # numeric gradient, so that
  gx <- (g(x + h, y) - g(x - h, y)) / (2 * h)   # swapping curves is a
  gy <- (g(x, y + h) - g(x, y - h)) / (2 * h)   # one-line change

  d <- abs(g(x, y)) / sqrt(gx^2 + gy^2)
  k <- which(is.finite(d) & d < tol)
  if (!length(k)) stop("no points within tol — loosen tol or widen the window")
  k <- k[sample(length(k), min(n, length(k)))]

  data.frame(x = x[k], y = y[k], d = d[k])
}

plot_variety <- function(pts, cols = c(gold, green, blue),
                         size = 0.5, alpha = 0.8, xlim, ylim) {
  ggplot(pts, aes(x, y, colour = d)) +
    geom_point(size = size, alpha = alpha) +
    scale_colour_gradientn(colours = cols) +
    coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) +
    bare
}

# some curves to try — the frame is 1.91:1, so wide curves fit best
lemniscate <- function(x, y, a = 3.2) (x^2 + y^2)^2 - a * (x^2 - y^2)
elliptic   <- function(x, y)          y^2 - x^3 + x        # sparse; lots of dead space
trifolium  <- function(x, y)          (x^2 + y^2)^2 - (x^3 - 3 * x * y^2)
astroid    <- function(x, y, a = 1)   (x^2 + y^2 - a^2)^3 + 27 * a^2 * x^2 * y^2


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

vx <- c(-2.5, 2.5); vy <- c(-1.31, 1.31)
save_og(
  plot_variety(sample_variety(lemniscate, c(-2.6, 2.6), c(-1.5, 1.5)),
               xlim = vx, ylim = vy),
  "og-a-variety.png"
)

save_og(plot_streams(field_spiral), "og-b-streams.png")

save_og(
  plot_hdr(rmix(), xlim = c(-3.6, 3.4), ylim = c(-1.7, 1.65)),
  "og-c-hdr.png"
)
