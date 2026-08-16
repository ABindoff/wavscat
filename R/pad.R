# Padding and unpadding.
#
# Scattering convolves circularly, because everything happens in the Fourier
# domain. Reflection padding to a power of two removes the resulting wrap-around
# artefacts and makes every transform length highly composite.

#' Split the required padding between the two ends of a signal
#'
#' @param N Padded length.
#' @param N_input Original length.
#' @return List with integer elements `left` and `right`.
#' @noRd
compute_padding <- function(N, N_input) {
  if (N < N_input) {
    stop("Padded length must be at least the signal length.", call. = FALSE)
  }
  to_add <- N - N_input
  pad_left <- to_add %/% 2L
  pad_right <- to_add - pad_left
  if (max(pad_left, pad_right) >= N_input) {
    stop(
      "Signal of length ", N_input, " is too short for the requested scale: ",
      "it would need ", max(pad_left, pad_right), " samples of reflection ",
      "padding, which is more than the signal itself. Reduce J or T.",
      call. = FALSE
    )
  }
  list(left = pad_left, right = pad_right)
}

#' Reflection padding
#'
#' Mirrors the signal at each end without repeating the edge sample, matching
#' `numpy.pad(mode = "reflect")`.
#'
#' @param x Numeric vector.
#' @param pad_left,pad_right Number of samples to add, each less than `length(x)`.
#' @return Numeric vector of length `length(x) + pad_left + pad_right`.
#' @noRd
pad_reflect <- function(x, pad_left, pad_right) {
  n <- length(x)
  if (pad_left >= n || pad_right >= n) {
    stop("Padding must be shorter than the signal.", call. = FALSE)
  }
  left <- if (pad_left > 0L) x[seq.int(pad_left + 1L, 2L)] else x[0L]
  right <- if (pad_right > 0L) x[seq.int(n - 1L, n - pad_right)] else x[0L]
  c(left, x, right)
}

#' Where the original signal sits inside the padded one, at every resolution
#'
#' After subsampling by `2^j` the boundaries move, and rounding outwards keeps
#' the retained window conservative.
#'
#' @param log2_T Log2 of the low-pass support.
#' @param J Log-scale of the transform.
#' @param i0 Zero-based start index of the signal within the padded signal.
#' @param i1 Zero-based exclusive end index.
#' @return List of integer vectors `start` and `end`, indexed by `j + 1`, both
#'   still zero-based.
#' @noRd
compute_border_indices <- function(log2_T, J, i0, i1) {
  m <- max(log2_T, J)
  start <- integer(m + 1L)
  end <- integer(m + 1L)
  start[1L] <- i0
  end[1L] <- i1
  for (j in seq_len(m)) {
    start[j + 1L] <- (start[j] %/% 2L) + (start[j] %% 2L)
    end[j + 1L] <- (end[j] %/% 2L) + (end[j] %% 2L)
  }
  list(start = start, end = end)
}

#' The indices of the original signal within a padded one, at level `j`
#'
#' @param borders Result of [compute_border_indices()].
#' @param j Subsampling level.
#' @return An integer vector of one-based indices.
#' @noRd
unpad_index <- function(borders, j) {
  seq.int(borders$start[j + 1L] + 1L, borders$end[j + 1L])
}

#' Cut a padded signal back to the extent of the original
#'
#' @param x Numeric vector at subsampling level `j`.
#' @param borders Result of [compute_border_indices()].
#' @param j Subsampling level.
#' @return Numeric vector.
#' @noRd
unpad <- function(x, borders, j) {
  x[unpad_index(borders, j)]
}
