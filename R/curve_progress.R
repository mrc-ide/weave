# =====================================================================
# curve_progress.R
# A compact, single-row progress bar drawn as a smooth sine wave that
# fills with a rainbow gradient as it grows.
#
# - One character row, drawn on a braille canvas (2x4 dots per cell),
#   so the wave reads as a continuous line.
# - The reached portion lights up red -> violet; the rest stays a
#   faint grey ghost so you can see the path ahead.
# - Stats (% and n/n) sit inline on the right, after the wave.
# - Redraw is a single carriage-return overwrite of one line, so there
#   is no cursor movement and therefore no flicker.
#
# Internal helper (not exported); used by gp_predict() for the
# posterior-draw loop. Requires a UTF-8 terminal with ANSI truecolor
# (any modern terminal, or the RStudio terminal pane).
#
# Usage:
#   pb <- make_curve_bar(total = 100)
#   for (i in 1:100) { Sys.sleep(0.02); pb$tick() }
#   pb$done()
# =====================================================================

# TRUE only in a terminal that can render the bar's ANSI truecolor + carriage-
# return redraw: an interactive session that is NOT a GUI console (Windows Rgui
# or macOS R.app emit raw escape codes and do not honour "\r" overwrites), and
# is either RStudio or a genuine tty.
ansi_tty <- function() {
  if (!interactive()) return(FALSE)
  if (.Platform$GUI %in% c("Rgui", "AQUA")) return(FALSE)
  Sys.getenv("RSTUDIO") == "1" || isatty(stdout())
}

make_curve_bar <- function(total,
                           width     = 48L,    # cells wide (the wave only)
                           cycles    = 10,     # number of pi-humps
                           thickness = 2L,     # line weight, in dots
                           hue_span  = 0.85,   # 0..1 of the colour wheel
                           stream    = stdout()) {
  width     <- as.integer(width)
  thickness <- max(1L, as.integer(thickness))
  H <- 4L                                       # one cell tall = 4 dot rows
  W <- width * 2L                               # dot columns

  # ---- Sine value v in [0,1] at full dot resolution -----------------
  x <- seq(0, cycles * pi, length.out = W)
  v <- (sin(x - pi / 2) + 1) / 2
  yf <- (1 - v) * (H - 1)                       # fractional row, 0 = top

  # ---- Rasterise a connected, thickened ribbon ----------------------
  grid <- matrix(FALSE, H, W)
  prev <- yf[1]
  for (xx in seq_len(W)) {
    cur <- yf[xx]
    lo  <- floor(min(prev, cur)); hi <- ceiling(max(prev, cur))
    extra <- thickness - (hi - lo + 1)
    if (extra > 0) { lo <- lo - extra %/% 2; hi <- hi + (extra - extra %/% 2) }
    lo <- max(0, lo); hi <- min(H - 1, hi)
    grid[(lo:hi) + 1L, xx] <- TRUE
    prev <- cur
  }

  # ---- Pack each cell's 2x4 dots into one braille glyph -------------
  leftbits  <- c(0x01, 0x02, 0x04, 0x40)
  rightbits <- c(0x08, 0x10, 0x20, 0x80)
  glyph <- character(width)
  for (cc in seq_len(width)) {
    dl <- (cc - 1L) * 2L + 1L; dr <- dl + 1L; bits <- 0L
    for (k in 1:4) {
      if (grid[k, dl]) bits <- bits + leftbits[k]
      if (dr <= W && grid[k, dr]) bits <- bits + rightbits[k]
    }
    glyph[cc] <- if (bits > 0L) intToUtf8(0x2800L + bits) else " "
  }

  # ---- Rainbow colour per cell --------------------------------------
  hsv2rgb <- function(h) {
    i <- floor(h * 6); f <- h * 6 - i; q <- 1 - f
    rgb <- switch(as.character(i %% 6),
                  "0" = c(1, f, 0), "1" = c(q, 1, 0), "2" = c(0, 1, f),
                  "3" = c(0, q, 1), "4" = c(f, 0, 1), "5" = c(1, 0, q),
                  c(1, 0, 0))
    round(rgb * 255)
  }
  hue     <- if (width > 1) seq(0, hue_span, length.out = width) else 0
  rgb     <- vapply(hue, hsv2rgb, numeric(3))
  rainbow <- sprintf("\033[1;38;2;%d;%d;%dm", rgb[1, ], rgb[2, ], rgb[3, ])
  ghost   <- "\033[38;2;78;78;88m"              # faint grey, unreached
  stat_c  <- "\033[1;97m"                        # bright white, the stats
  reset   <- "\033[0m"

  # Fixed-width stats (constant length => nothing jitters): "  50%  70/140"
  total   <- as.integer(total)
  digits  <- nchar(as.character(total))
  statfmt <- sprintf("  %%3d%%%%  %%%dd/%%d", digits)

  st <- new.env(parent = emptyenv())
  st$current <- 0L; st$total <- total

  render <- function() {
    p        <- max(0, min(1, st$current / st$total))
    revealed <- floor(p * width)
    cells <- vapply(seq_len(width), function(cc)
      paste0(if (cc <= revealed) rainbow[cc] else ghost, glyph[cc]),
      character(1))
    stats <- sprintf(statfmt, round(p * 100), st$current, total)
    cat("\r", paste0(cells, collapse = ""), reset, stat_c, stats,
        reset, "\033[K", sep = "", file = stream)
    flush(stream)
  }

  list(
    tick = function(by = 1L) {
      st$current <- min(st$total, st$current + as.integer(by)); render(); invisible()
    },
    set = function(value) {
      st$current <- max(0L, min(st$total, as.integer(value))); render(); invisible()
    },
    done = function() {
      st$current <- st$total; render(); cat("\n", file = stream); flush(stream); invisible()
    }
  )
}
