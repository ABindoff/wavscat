# Array primitives for the scattering cascade.
#
# Every elementwise / FFT / resampling operation the cascade performs goes
# through this file. That is deliberate: it is the single seam at which a
# compiled engine (extendr + rustfft) can be substituted without touching the
# filter bank, the path enumeration or the user-facing API. Keep these
# functions free of scattering-specific logic.

#' Forward discrete Fourier transform
#'
#' Unnormalised, matching the convention of `numpy.fft.fft`.
#'
#' @param x Numeric or complex vector.
#' @return Complex vector of the same length.
#' @noRd
bk_fft <- function(x) {
  stats::fft(x)
}

#' Inverse discrete Fourier transform
#'
#' Normalised by `1 / length(x)`, matching `numpy.fft.ifft`. Base R's
#' `fft(inverse = TRUE)` omits that factor.
#'
#' @param x Complex vector.
#' @return Complex vector of the same length.
#' @noRd
bk_ifft <- function(x) {
  stats::fft(x, inverse = TRUE) / length(x)
}

#' Inverse transform of a spectrum known to represent a real signal
#'
#' @param x Complex vector.
#' @return Numeric vector.
#' @noRd
bk_irfft <- function(x) {
  Re(bk_ifft(x))
}

#' Periodise a spectrum, which subsamples the corresponding signal
#'
#' Subsampling a signal by `k` in time is periodisation by `k` in frequency.
#' Element `n` of the result reduces `x[i * m + n]` over `i = 0, ..., k - 1`,
#' where `m = length(x) / k`.
#'
#' The reduction differs by use. Subsampling an actual signal spectrum averages,
#' which carries the `1 / k` that keeps the inverse transform correctly scaled.
#' Rescaling a *filter* to a coarser resolution sums, so that its response
#' against an already-subsampled signal matches its response at full resolution.
#'
#' @param x Numeric or complex vector whose length is divisible by `k`.
#' @param k Positive integer subsampling factor.
#' @param reduce Either `"mean"` (signals) or `"sum"` (filters).
#' @return Vector of length `length(x) / k`.
#' @noRd
bk_periodize <- function(x, k, reduce = c("mean", "sum")) {
  k <- as.integer(k)
  if (k <= 1L) {
    return(x)
  }
  n <- length(x)
  if (n %% k != 0L) {
    stop("Cannot periodise a vector of length ", n, " by ", k, ".", call. = FALSE)
  }
  reduce <- match.arg(reduce)
  m <- n %/% k

  # Halving is the common case on the longest vectors, and a single vectorised
  # add beats any matrix reduction there by a wide margin.
  if (k == 2L) {
    acc <- x[seq_len(m)] + x[m + seq_len(m)]
    return(if (reduce == "mean") acc / 2 else acc)
  }

  # Column-major fill puts x[i * m + n] at position [n, i], so folding over
  # blocks is a row reduction. Expressing it as a matrix-vector product hands
  # the work to BLAS, which is faster than rowMeans on complex input.
  dim(x) <- c(m, k)
  weight <- if (reduce == "mean") rep(1 / k, k) else rep(1, k)
  as.vector(x %*% weight)
}

#' Complex modulus
#' @noRd
bk_modulus <- function(x) {
  Mod(x)
}

# ---------------------------------------------------------------------------
# Matrix primitives, used by joint time-frequency scattering.
#
# Joint scattering carries a matrix indexed by [frequency band, time] and
# convolves along each axis in turn. Base R's mvfft transforms the columns of a
# matrix in one C-level call, so the frequency direction is free; the time
# direction costs two transposes.
# ---------------------------------------------------------------------------

#' Forward transform down each column, that is along the frequency axis
#' @noRd
bk_fft_freq <- function(M) {
  stats::mvfft(M)
}

#' @noRd
bk_ifft_freq <- function(M) {
  stats::mvfft(M, inverse = TRUE) / nrow(M)
}

#' Periodise a matrix along one axis
#'
#' The matrix analogue of [bk_periodize()]: subsamples the frequency axis when
#' `axis` is 1 and the time axis when it is 2.
#'
#' @param M A matrix.
#' @param k Subsampling factor; must divide the length of the chosen axis.
#' @param axis 1 for rows (frequency), 2 for columns (time).
#' @param reduce `"mean"` for signals, `"sum"` for filters.
#' @noRd
bk_periodize_axis <- function(M, k, axis, reduce = c("mean", "sum")) {
  k <- as.integer(k)
  if (k <= 1L) {
    return(M)
  }
  reduce <- match.arg(reduce)
  nr <- nrow(M)
  nc <- ncol(M)

  # Halving avoids building the intermediate three-way array at all.
  if (k == 2L) {
    out <- if (axis == 1L) {
      M[seq_len(nr %/% 2L), , drop = FALSE] +
        M[nr %/% 2L + seq_len(nr %/% 2L), , drop = FALSE]
    } else {
      M[, seq_len(nc %/% 2L), drop = FALSE] +
        M[, nc %/% 2L + seq_len(nc %/% 2L), drop = FALSE]
    }
    return(if (reduce == "mean") out / 2 else out)
  }

  # Both cases are a pure reshape: column-major storage already lays the blocks
  # out contiguously along the axis being folded.
  if (axis == 1L) {
    if (nr %% k != 0L) {
      stop("Cannot periodise ", nr, " rows by ", k, ".", call. = FALSE)
    }
    a <- array(M, dim = c(nr %/% k, k, nc))
    acc <- a[, 1L, , drop = FALSE]
    for (i in seq_len(k - 1L) + 1L) acc <- acc + a[, i, , drop = FALSE]
    out <- matrix(acc, nrow = nr %/% k, ncol = nc)
  } else {
    if (nc %% k != 0L) {
      stop("Cannot periodise ", nc, " columns by ", k, ".", call. = FALSE)
    }
    a <- array(M, dim = c(nr, nc %/% k, k))
    acc <- a[, , 1L, drop = FALSE]
    for (i in seq_len(k - 1L) + 1L) acc <- acc + a[, , i, drop = FALSE]
    out <- matrix(acc, nrow = nr, ncol = nc %/% k)
  }
  if (reduce == "mean") out / k else out
}

#' Pointwise product of a spectrum with a filter
#' @noRd
bk_cdgmm <- function(a, b) {
  a * b
}

#' Sum over time, used when averaging is global rather than local
#' @noRd
bk_average_global <- function(x) {
  sum(x)
}
