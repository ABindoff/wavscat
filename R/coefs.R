# Path metadata and the printing of operators and coefficients.

#' Describe the coefficients an operator produces
#'
#' One row per scattering path, in the order the coefficients are returned.
#' Frequencies are reported both in cycles per sample, where 0.5 is Nyquist, and
#' in Hz when the operator knows its sampling rate.
#'
#' The second-order centre frequency `xi2` is a *modulation* frequency: it is
#' the rate at which the envelope of the first-order band `xi1` fluctuates. A
#' path with `xi1` near 1 kHz and `xi2` near 5 Hz describes energy around 1 kHz
#' being amplitude-modulated five times a second, which is roughly the syllable
#' rate of speech.
#'
#' @param op A scattering operator, or a `wavscat_coefs` object.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{path}{Label identifying the coefficient.}
#'     \item{order}{Scattering order: 0, 1 or 2.}
#'     \item{n1, n2}{Indices of the wavelets used, `NA` where not applicable.}
#'     \item{j1, j2}{Dyadic scales of those wavelets.}
#'     \item{xi1, xi2}{Centre frequencies in cycles per sample.}
#'     \item{sigma1, sigma2}{Bandwidths in cycles per sample.}
#'     \item{xi1_hz, xi2_hz}{Centre frequencies in Hz, present only when `sr`
#'       was supplied.}
#'   }
#'
#' @examples
#' op <- scattering_1d(n = 4096, J = 6, Q = 8, sr = 16000)
#' m <- scat_meta(op)
#' head(m)
#' table(m$order)
#' @export
scat_meta <- function(op) {
  if (inherits(op, "wavscat_coefs")) {
    return(op$meta)
  }
  if (!inherits(op, "wavscat_op")) {
    stop("op must be a scattering operator or a wavscat_coefs object.",
         call. = FALSE)
  }
  paths <- op$paths
  # Put the identifying columns first, then whatever the operator produced.
  lead <- c("path", "order", "n1", "n2", "n_fr")
  cols <- c(intersect(lead, names(paths)), setdiff(names(paths), lead))
  out <- tibble::as_tibble(paths[, cols, drop = FALSE])

  if (!is.null(op$sr)) {
    out$xi1_hz <- out$xi1 * op$sr
    out$xi2_hz <- out$xi2 * op$sr
  }
  if ("xi_fr" %in% names(out)) {
    # xi_fr counts cycles per first-order band; there are Q1 bands per octave,
    # so this is the more interpretable rate of spectral drift.
    out$xi_fr_octave <- out$xi_fr * op$Q[1L]
  }
  out
}

#' Number of coefficients an operator produces
#'
#' @param op A scattering operator.
#' @param by_order Whether to break the count down by scattering order.
#' @return A single integer, or a named integer vector when `by_order` is `TRUE`.
#' @examples
#' op <- scattering_1d(n = 4096, J = 6, Q = 8)
#' scat_output_size(op)
#' scat_output_size(op, by_order = TRUE)
#' @export
scat_output_size <- function(op, by_order = FALSE) {
  if (!inherits(op, "wavscat_op")) {
    stop("op must be a scattering operator.", call. = FALSE)
  }
  if (!by_order) {
    return(nrow(op$paths))
  }
  counts <- table(factor(op$paths$order, levels = 0:op$max_order))
  stats::setNames(as.integer(counts), paste0("order", names(counts)))
}

#' @export
print.wavscat_op <- function(x, ...) {
  cat("<wavscat_op>",
      if (x$type == "time") {
        "time scattering"
      } else {
        sprintf("joint time-frequency scattering (format = \"%s\")", x$format)
      }, "\n")
  cat("  signal length : ", x$n,
      if (!is.null(x$sr)) sprintf(" samples (%.4g s at %g Hz)",
                                  x$n / x$sr, x$sr) else " samples", "\n", sep = "")
  cat("  J = ", x$J, "   Q = ", paste(x$Q, collapse = ", "),
      "   max order = ", x$max_order, "\n", sep = "")
  if (identical(x$type, "jtfs")) {
    cat("  J_fr = ", x$J_fr, "   Q_fr = ", x$Q_fr, "   F = ",
        if (identical(x$average_fr, "global")) "global"
        else if (isFALSE(x$average_fr)) "none" else x$F,
        " bands\n", sep = "")
  }

  avg <- if (identical(x$average, "global")) {
    "global (one coefficient per path)"
  } else if (isFALSE(x$average)) {
    "none"
  } else if (!is.null(x$sr)) {
    sprintf("T = %g samples (%.4g s)", x$T, x$T / x$sr)
  } else {
    sprintf("T = %g samples", x$T)
  }
  cat("  averaging     : ", avg, "\n", sep = "")

  if (identical(x$average, "local")) {
    n_time <- x$borders$end[x$log2_stride + 1L] -
      x$borders$start[x$log2_stride + 1L]
    sr_out <- out_sample_rate(x)
    cat("  output        : ", nrow(x$paths), " paths x ", n_time, " time points",
        if (!is.na(sr_out)) sprintf(" (%.4g Hz)", sr_out) else "", "\n", sep = "")
  } else {
    cat("  output        : ", nrow(x$paths), " paths\n", sep = "")
  }

  counts <- scat_output_size(x, by_order = TRUE)
  cat("  paths by order: ", paste(names(counts), counts, sep = " = ",
                                  collapse = ", "), "\n", sep = "")
  if (!is.null(x$sr)) {
    xi1 <- x$paths$xi1[x$paths$order == 1L]
    cat(sprintf("  band-pass     : %.4g to %.4g Hz\n",
                min(xi1) * x$sr, max(xi1) * x$sr))
  }
  invisible(x)
}

#' @export
print.wavscat_coefs <- function(x, ...) {
  cat("<wavscat_coefs>\n")
  if (is.array(x$coef)) {
    d <- dim(x$coef)
    if (length(d) == 4L) {
      cat("  ", d[1L], " paths x ", d[2L], " bands x ", d[3L],
          " time points x ", d[4L], " channel",
          if (d[4L] != 1L) "s" else "", "\n", sep = "")
    } else {
      cat("  ", d[1L], " paths x ", d[2L], " time points x ", d[3L],
          " channel", if (d[3L] != 1L) "s" else "", "\n", sep = "")
    }
  } else {
    cat("  ", length(x$coef), " paths (list form)\n", sep = "")
  }
  if (!is.na(x$spec$sr_out)) {
    cat(sprintf("  output rate   : %.4g Hz\n", x$spec$sr_out))
  }
  cat("  paths by order: ",
      paste(names(table(x$meta$order)), as.integer(table(x$meta$order)),
            sep = " = ", collapse = ", "), "\n", sep = "")
  if (length(x$transforms)) {
    cat("  applied       : ", paste(x$transforms, collapse = ", "), "\n", sep = "")
  }
  if (is.array(x$coef)) {
    e <- scat_energy(x)
    cat(sprintf("  energy by order: %s\n",
                paste(sprintf("%d = %.3g", e$order, e$energy), collapse = ", ")))
  }
  invisible(x)
}

#' Energy carried by each scattering order
#'
#' The squared L2 norm of the coefficients, summed within each order and
#' averaged over channels. Useful for checking that a choice of `J`, `Q` and `T`
#' actually captures the structure of a signal: if order 2 holds almost no
#' energy, its paths are not earning their place in the feature vector.
#'
#' @param x A `wavscat_coefs` object with array storage.
#' @return A tibble with columns `order` and `energy`.
#' @examples
#' S <- scattering(as.numeric(arima.sim(list(ar = 0.9), 2048)), J = 6, Q = 8)
#' scat_energy(S)
#' @export
scat_energy <- function(x) {
  if (!inherits(x, "wavscat_coefs")) {
    stop("x must be a wavscat_coefs object.", call. = FALSE)
  }
  if (!is.array(x$coef)) {
    stop("scat_energy() needs array storage; rebuild with out_type = \"array\".",
         call. = FALSE)
  }
  n_channels <- utils::tail(dim(x$coef), 1L)
  e2 <- apply(x$coef^2, 1L, sum) / n_channels
  agg <- tapply(e2, x$meta$order, sum)
  tibble::tibble(
    order = as.integer(names(agg)),
    energy = as.numeric(agg)
  )
}
