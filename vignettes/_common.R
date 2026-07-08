# =============================================================================
# Shared styling for the weave vignettes.
#
# Sourced once at the top of each vignette (the working directory during the
# vignette build is `vignettes/`, so `source("_common.R")` is safe). It sets the
# knitr chunk defaults and defines a single ggplot2 theme + palette so every
# figure across the vignettes shares one visual identity.
#
# The palette is drawn from the hummingbird photograph used in the gentle
# introduction (man/figures/Hummingbird.jpg): pine, teal, rose, moss, cream.
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
# The palette is drawn from the hummingbird photograph that opens the gentle
# introduction (man/figures/Hummingbird.jpg): deep pine plumage, an iridescent
# teal back, hibiscus rose, moss-green leaves, and a cream bokeh background.
weave_navy   <- "#2c463f"  # deep pine: headings / axis text / the truth
weave_linen  <- "#f7f4ea"  # cream bokeh, for panels and strips
weave_grid   <- "#e4e7d9"  # faint sage grid line

# Categorical palette for GP draws etc., in the same photograph's hues.
weave_pal <- c(
  teal  = "#2f8f77",  # iridescent back
  rose  = "#d4699e",  # hibiscus petals
  moss  = "#7ba05b",  # leaves
  honey = "#c9973b",  # warm bokeh light
  pine  = "#2c463f",
  sage  = "#a7bfae"
)

# Stable, named colours so the same role keeps the same colour across plots.
weave_cols <- c(
  truth      = weave_navy,
  prediction = "#2f8f77",   # teal
  estimate   = "#2f8f77",
  held_out   = "#d4699e",   # rose
  observed   = "#8f8779",   # warm taupe (the bird's breast)
  kernel     = "#b23a6f"    # deep rose (the flower's centre)
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
