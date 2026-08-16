# Diagnostic plots.
#
# Both are for deciding whether a choice of J, Q and T suits a signal, so they
# favour legibility over decoration: recessive axes, a single-hue ramp ordered
# by frequency rather than arbitrary categorical colours, and a log frequency
# axis, which is the geometry a constant-Q bank actually lives in.

#' Plot the frequency response of a filter bank
#'
#' Shows every wavelet in the bank against frequency, with the averaging filter
#' drawn separately. Read it to check two things: that the wavelets tile the
#' band you care about without gaps, and that the low-pass sits where you expect
#' given your choice of `T`.
#'
#' Colour runs light to dark with centre frequency; it encodes position in the
#' bank, not identity, so individual curves are deliberately not distinguishable.
#'
#' @param x A scattering operator.
#' @param order Which filter bank to draw, 1 or 2.
#' @param lp Whether to overlay the Littlewood-Paley sum, that is the total
#'   power the bank applies at each frequency. A flat sum near one means the bank
#'   covers the axis evenly; dips mean frequencies that the transform sees only
#'   faintly.
#' @param n_freq Number of frequency points to draw.
#' @param col_ramp Name of an [grDevices::hcl.colors()] palette for the wavelets.
#' @param ... Passed to [graphics::plot()].
#'
#' @return `x`, invisibly.
#'
#' @examples
#' op <- scattering_1d(n = 8192, J = 8, Q = 8, sr = 16000)
#' plot(op)
#' plot(op, lp = TRUE)
#' @export
plot.wavscat_op <- function(x, order = 1, lp = FALSE, n_freq = 2048,
                            col_ramp = "Teal", ...) {
  order <- as.integer(order)
  if (!order %in% c(1L, 2L)) {
    stop("order must be 1 or 2.", call. = FALSE)
  }
  bank <- if (order == 1L) x$filters$psi1 else x$filters$psi2

  N <- x$n_padded
  freq <- fftfreq(N)
  keep <- which(freq > 0)
  # Thin to a manageable number of points, keeping the log spacing sensible.
  if (length(keep) > n_freq) {
    keep <- keep[unique(round(exp(seq(0, log(length(keep)), length.out = n_freq))))]
    keep <- keep[keep >= 1]
  }
  f <- freq[keep]

  in_hz <- !is.null(x$sr)
  fx <- if (in_hz) f * x$sr else f
  xlab <- if (in_hz) "frequency (Hz)" else "frequency (cycles per sample)"

  curves <- lapply(bank, function(p) p$levels[[1L]][keep])
  phi <- x$filters$phi$levels[[1L]][keep]
  ymax <- max(1, vapply(curves, max, 0))

  ink <- "grey20"
  muted <- "grey45"
  grid_col <- "grey88"

  op_par <- graphics::par(mar = c(4.1, 4.1, 3.1, 1.1))
  on.exit(graphics::par(op_par), add = TRUE)

  graphics::plot(
    range(fx), c(0, ymax * 1.02), type = "n", log = "x",
    xlab = xlab, ylab = "gain", axes = FALSE,
    main = sprintf("wavscat filter bank, order %d  (J = %d, Q = %d)",
                   order, x$J, x$Q[order]),
    col.main = ink, col.lab = muted, ...
  )
  graphics::grid(col = grid_col, lty = 1, lwd = 0.5)
  graphics::axis(1, col = grid_col, col.axis = muted, lwd = 0.5)
  graphics::axis(2, col = grid_col, col.axis = muted, lwd = 0.5, las = 1)

  # Light to dark with increasing centre frequency: a magnitude encoding, not an
  # identity one.
  cols <- grDevices::hcl.colors(length(curves), col_ramp, rev = TRUE)
  ord <- order(vapply(bank, `[[`, 0, "xi"))
  for (i in seq_along(ord)) {
    graphics::lines(fx, curves[[ord[i]]], col = cols[i], lwd = 1.1)
  }
  graphics::lines(fx, phi, col = ink, lwd = 2)

  if (lp) {
    lpsum <- phi^2
    for (h in curves) lpsum <- lpsum + h^2
    graphics::lines(fx, lpsum, col = "#B4451F", lwd = 2, lty = 2)
  }

  graphics::legend(
    "topleft", bty = "n", cex = 0.85, text.col = muted,
    legend = c(sprintf("%d wavelets", length(curves)),
               sprintf("low-pass (T = %g)", x$T),
               if (lp) "Littlewood-Paley sum"),
    col = c(cols[length(cols) %/% 2L], ink, if (lp) "#B4451F"),
    lwd = c(1.1, 2, if (lp) 2), lty = c(1, 1, if (lp) 2)
  )
  invisible(x)
}


#' Plot scattering coefficients as a scalogram
#'
#' Displays one order of coefficients as an image: time across, frequency band up,
#' magnitude as colour. For order 1 this is a scalogram, directly comparable to a
#' constant-Q spectrogram. For order 2 the rows are `(n1, n2)` pairs grouped by
#' their first-order parent, with faint rules marking the group boundaries.
#'
#' Log-compressed coefficients are usually far more readable than raw ones; see
#' [scat_log()].
#'
#' @param x A `wavscat_coefs` object with array storage.
#' @param order Scattering order to display.
#' @param channel Channel index or name.
#' @param col_ramp Name of an [grDevices::hcl.colors()] palette. The default is
#'   perceptually uniform and safe for colour vision deficiency.
#' @param ... Passed to [graphics::image()].
#'
#' @return `x`, invisibly.
#'
#' @examples
#' t <- seq_len(8192)
#' x <- sin(2 * pi * (0.002 + 0.2 * t / 8192) * t)
#' S <- scattering(x, J = 8, Q = 8, sr = 8000)
#' plot(scat_log(S, eps = scat_eps(S)))
#' @export
plot.wavscat_coefs <- function(x, order = 1, channel = 1, col_ramp = "viridis",
                               ...) {
  check_coefs(x)
  if (!is.array(x$coef)) {
    stop("Plotting needs array storage; rebuild the operator with ",
         "out_type = \"array\".", call. = FALSE)
  }
  if (length(dim(x$coef)) == 4L) {
    stop("Joint coefficients are one image per path, so there is no single ",
         "scalogram to draw. Build the operator with format = \"time\", or ",
         "image one path yourself from the [path, band, time, channel] array.",
         call. = FALSE)
  }
  ch <- if (is.character(channel)) match(channel, x$channels) else as.integer(channel)
  if (is.na(ch) || ch < 1L || ch > length(x$channels)) {
    stop("channel must name or index one of: ",
         paste(x$channels, collapse = ", "), ".", call. = FALSE)
  }

  rows <- which(x$meta$order == order)
  if (length(rows) == 0L) {
    stop("This object holds no order-", order, " coefficients.", call. = FALSE)
  }
  if (dim(x$coef)[2L] < 2L) {
    stop("There is only one time point per path, so there is nothing to plot ",
         "as a scalogram. Rebuild without T = \"global\".", call. = FALSE)
  }

  # Order rows by band so that frequency increases up the image.
  rows <- rows[order(x$meta$xi1[rows], x$meta$xi2[rows], na.last = FALSE)]
  m <- x$coef[rows, , ch]

  n_time <- ncol(m)
  sr_out <- x$spec$sr_out
  tvals <- if (!is.na(sr_out)) (seq_len(n_time) - 1) / sr_out else seq_len(n_time)
  xlab <- if (!is.na(sr_out)) "time (s)" else "time (samples of the output)"

  ink <- "grey20"
  muted <- "grey45"

  op_par <- graphics::par(mar = c(4.1, 4.6, 3.1, 1.1))
  on.exit(graphics::par(op_par), add = TRUE)

  graphics::image(
    x = tvals, y = seq_along(rows), z = t(m),
    col = grDevices::hcl.colors(256, col_ramp),
    xlab = xlab, ylab = "", axes = FALSE,
    main = sprintf("order %d%s", order,
                   if (length(x$transforms)) {
                     paste0("  (", paste(x$transforms, collapse = ", "), ")")
                   } else ""),
    col.main = ink, col.lab = muted, ...
  )
  graphics::axis(1, col = "grey70", col.axis = muted, lwd = 0.5)

  if (order == 1L) {
    # Label a handful of bands rather than all of them.
    at <- unique(round(seq(1, length(rows), length.out = 6)))
    lab <- if (!is.null(x$spec$sr)) {
      sprintf("%.3g", x$meta$xi1[rows][at] * x$spec$sr)
    } else {
      sprintf("%.3g", x$meta$xi1[rows][at])
    }
    graphics::axis(2, at = at, labels = lab, las = 1, col = "grey70",
                   col.axis = muted, lwd = 0.5)
    graphics::mtext(if (!is.null(x$spec$sr)) "band (Hz)" else "band (cycles/sample)",
                    side = 2, line = 3.2, col = muted, cex = 0.9)
  } else {
    # Rule off each first-order parent group.
    grp <- x$meta$n1[rows]
    bnd <- which(diff(grp) != 0) + 0.5
    graphics::abline(h = bnd, col = grDevices::adjustcolor("white", 0.25), lwd = 0.5)
    graphics::mtext("(band, modulation) pair", side = 2, line = 1.4,
                    col = muted, cex = 0.9)
  }
  graphics::box(col = "grey70", lwd = 0.5)
  invisible(x)
}
