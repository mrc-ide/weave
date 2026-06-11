# =============================================================================
# Shared styling for the weave vignettes.
#
# Sourced once at the top of each vignette (the working directory during the
# vignette build is `vignettes/`, so `source("_common.R")` is safe). It sets the
# knitr chunk defaults and defines a single ggplot2 theme + palette so every
# figure across the vignettes shares one visual identity.
#
# The palette is taken from the hex logo (man/figures/Weave.png): a dark-navy
# border, a warm woven-linen background, and the vivid multicolour lettering.
# =============================================================================

knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 7,
  fig.height = 4.2,
  dpi    = 150,
  dev    = "png",
  fig.align = "center",
  warning = FALSE,
  message = FALSE
)

# --- brand colours -----------------------------------------------------------
weave_navy   <- "#1b2a4a"  # logo border / headings / axis text
weave_linen  <- "#f7f1e8"  # logo background, softened for panels and strips
weave_grid   <- "#e7ded0"  # faint warm grid line

# Vivid lettering, used as the categorical palette for GP draws etc.
weave_pal <- c(
  violet  = "#7b2ff7",
  magenta = "#ff2d95",
  blue    = "#2e7fff",
  green   = "#1dd1a1",
  gold    = "#feca57",
  sky     = "#48dbfb"
)

# Stable, named colours so the same role keeps the same colour across plots.
weave_cols <- c(
  truth      = weave_navy,
  prediction = "#2e7fff",   # logo blue
  estimate   = "#2e7fff",
  held_out   = "#ff2d95",   # logo magenta
  observed   = "grey55",
  kernel     = "#7b2ff7"    # logo violet
)

# --- one theme for every figure ----------------------------------------------
theme_weave <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text             = ggplot2::element_text(colour = weave_navy),
      plot.title       = ggplot2::element_text(face = "bold", size = base_size + 3,
                                               colour = weave_navy),
      plot.subtitle    = ggplot2::element_text(colour = "grey35", size = base_size - 1),
      axis.title       = ggplot2::element_text(face = "bold", colour = weave_navy),
      axis.text        = ggplot2::element_text(colour = "grey30"),
      panel.grid.major = ggplot2::element_line(colour = weave_grid, linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line        = ggplot2::element_line(colour = weave_navy, linewidth = 0.4),
      legend.title     = ggplot2::element_text(face = "bold", colour = weave_navy),
      legend.text      = ggplot2::element_text(colour = weave_navy),
      legend.position  = "top",
      strip.background = ggplot2::element_rect(fill = weave_linen, colour = NA),
      strip.text       = ggplot2::element_text(face = "bold", colour = weave_navy,
                                               size = base_size - 3),
      plot.margin      = ggplot2::margin(8, 12, 8, 8)
    )
}

# Convenience scales tied to the brand palette.
scale_colour_weave <- function(...) {
  ggplot2::scale_colour_manual(values = unname(weave_pal), ...)
}
scale_fill_weave <- function(...) {
  ggplot2::scale_fill_manual(values = unname(weave_pal), ...)
}
