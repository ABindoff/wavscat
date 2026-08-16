# Construction of a one-dimensional time scattering operator.

#' Create a 1D wavelet scattering operator
#'
#' Builds the Morlet filter banks and the path table for a time scattering
#' transform of signals of a fixed length. Building the operator is the
#' expensive step; apply it to as many signals as you like with
#' [scat_transform()].
#'
#' @section Choosing the parameters:
#'
#' `J` sets the largest wavelet scale, `2^J` samples, and so the lowest
#' frequency the bank reaches, roughly `sr / 2^J` Hz. `Q` sets the frequency
#' resolution: `Q` wavelets per octave, each with quality factor about `Q`.
#' `Q = 8` at the first order resolves close-together frequencies and suits
#' audio; `Q = 1` gives wide filters that localise transients, which is the
#' usual second-order choice.
#'
#' `T` is separate from `J` and is the one that most often needs attention. It
#' is the support of the averaging filter, in samples, and therefore sets both
#' the scale of translation invariance and the output sampling interval. It
#' defaults to `2^J`, which for short recordings can average away everything you
#' care about: an 1800-sample tapping trace at 30 Hz with `J = 8` would be
#' averaged over `2^8 / 30`, about 8.5 seconds. Use `T_sec` to set it in
#' seconds instead, or `T = "global"` to reduce each path to a single number
#' per signal.
#'
#' @param n Length of the signals to be transformed, in samples.
#' @param J Log-scale of the transform. The largest wavelet has support `2^J`
#'   samples, so `J` must satisfy `2^J <= n`.
#' @param Q Wavelets per octave. Either one number, applied to the first order
#'   with the second order set to 1, or two numbers `c(Q1, Q2)`.
#' @param T Support of the low-pass averaging filter in samples. `NULL`, the
#'   default, uses `2^J`. `"global"` averages each path over the whole signal,
#'   giving one coefficient per path. `0` disables averaging, in which case
#'   `out_type` must be `"list"`.
#'
#'   Be careful with `"global"`. It sums over the reflection-padded signal, so
#'   samples near each end are counted twice, matching the behaviour of the
#'   Kymatio reference implementation. To get one value per path without that
#'   distortion, use a local `T` and reduce the time axis afterwards with
#'   `scat_features(summary = "mean")`.
#' @param T_sec Averaging support in seconds. Requires `sr`, and overrides `T`.
#' @param max_order Highest scattering order, 1 or 2. Order 2 is where most of
#'   the discriminative power beyond a mel spectrogram lies.
#' @param stride Output subsampling factor, a power of two no greater than `T`.
#'   Defaults to `T` rounded down to a power of two, which is critical sampling.
#'   A smaller stride oversamples the output in time.
#' @param sr Sampling rate in Hz. Optional, and used only to report centre
#'   frequencies in Hz and to interpret `T_sec`; it does not change the
#'   coefficients.
#' @param out_type `"array"` to stack paths into a matrix, or `"list"` to keep
#'   each path separate. Only `"list"` is possible when `T = 0`.
#'
#' @return An object of class `wavscat_op`.
#'
#' @references
#' Anden, J. and Mallat, S. (2014). Deep scattering spectrum.
#' *IEEE Transactions on Signal Processing*, 62(16), 4114-4128.
#' \doi{10.1109/TSP.2014.2326991}
#'
#' @examples
#' op <- scattering_1d(n = 4096, J = 6, Q = 8)
#' op
#'
#' # A tapping trace at 30 Hz, averaged over one second rather than over 2^J
#' # samples, which here would be more than eight seconds.
#' op_tap <- scattering_1d(n = 1800, J = 8, Q = c(8, 1), T_sec = 1, sr = 30)
#' op_tap
#' @export
scattering_1d <- function(n, J, Q = 8, T = NULL, T_sec = NULL, max_order = 2,
                          stride = NULL, sr = NULL,
                          out_type = c("array", "list")) {
  out_type <- match.arg(out_type)

  n <- check_count(n, "n")
  J <- check_count(J, "J")
  if (2^J > n) {
    stop(
      "J = ", J, " implies a maximum wavelet support of 2^", J, " = ", 2^J,
      " samples, which exceeds the signal length n = ", n, ". Use J <= ",
      floor(log2(n)), ".",
      call. = FALSE
    )
  }
  Q <- parse_Q(Q)
  max_order <- as.integer(max_order)
  if (!max_order %in% c(1L, 2L)) {
    stop("max_order must be 1 or 2, got ", max_order, ".", call. = FALSE)
  }
  if (!is.null(sr) && (!is.numeric(sr) || length(sr) != 1L || sr <= 0)) {
    stop("sr must be a single positive number, or NULL.", call. = FALSE)
  }

  if (!is.null(T_sec)) {
    if (is.null(sr)) {
      stop("T_sec needs sr, so that seconds can be converted to samples.",
           call. = FALSE)
    }
    T <- round(T_sec * sr)
  }

  parsed <- parse_T(T, J, n)
  T <- parsed$T
  average <- parsed$average
  log2_T <- as.integer(floor(log2(T)))

  if (average == FALSE && out_type == "array") {
    stop("With T = 0 the paths have different lengths, so out_type must be ",
         "\"list\".", call. = FALSE)
  }

  log2_stride <- parse_stride(stride, log2_T, average)

  r_psi <- sqrt(0.5)
  sigma0 <- 0.1
  alpha <- 5

  # Pad by three times the half support of the averaging filter, which is where
  # its tail has fallen below one part in a thousand, but never so far that the
  # reflection would have to repeat.
  min_to_pad <- 3 * compute_temporal_support(gauss_1d(n, sigma0 / T))
  J_max_support <- floor(log2(3 * n - 2))
  J_pad <- min(ceiling(log2(n + 2 * min_to_pad)), J_max_support)
  n_padded <- as.integer(2^J_pad)

  padding <- compute_padding(n_padded, n)
  borders <- compute_border_indices(
    log2_T, J, padding$left, padding$left + n
  )

  filters <- scattering_filter_factory(n_padded, J, Q, T, alpha, r_psi, sigma0)

  op <- structure(
    list(
      type = "time",
      n = n, J = J, Q = Q, T = T, max_order = max_order,
      average = average, log2_T = log2_T, log2_stride = log2_stride,
      sr = sr, out_type = out_type,
      n_padded = n_padded, pad_left = padding$left, pad_right = padding$right,
      borders = borders, filters = filters,
      r_psi = r_psi, sigma0 = sigma0, alpha = alpha
    ),
    class = "wavscat_op"
  )
  op$paths <- enumerate_paths(op)
  op
}

#' Wavelets per octave, normalised to a pair
#' @noRd
parse_Q <- function(Q) {
  if (!is.numeric(Q) || length(Q) < 1L || length(Q) > 2L) {
    stop("Q must be one or two numbers.", call. = FALSE)
  }
  if (any(Q < 1)) {
    stop("Q must be at least 1, got ", paste(Q, collapse = ", "), ".",
         call. = FALSE)
  }
  if (any(Q != round(Q))) {
    stop("Q must be whole numbers.", call. = FALSE)
  }
  Q <- as.integer(Q)
  if (length(Q) == 1L) c(Q, 1L) else Q
}

#' Averaging support and averaging mode
#'
#' @return List with `T` in samples and `average` one of `"local"`, `"global"`
#'   or `FALSE`.
#' @noRd
parse_T <- function(T, J, n, arg = "T") {
  if (is.null(T)) {
    return(list(T = 2^J, average = "local"))
  }
  if (identical(T, "global")) {
    return(list(T = 2^J, average = "global"))
  }
  if (!is.numeric(T) || length(T) != 1L) {
    stop(arg, " must be a single number, 0, \"global\", or NULL.", call. = FALSE)
  }
  if (T == 0) {
    return(list(T = 2^J, average = FALSE))
  }
  if (T < 1) {
    stop(arg, " must be 0 or at least 1, got ", T, ".", call. = FALSE)
  }
  if (T > n) {
    stop(
      arg, " = ", T, " exceeds the signal length ", n,
      ". For averaging over the whole signal use ", arg, " = \"global\".",
      call. = FALSE
    )
  }
  list(T = T, average = "local")
}

#' @noRd
parse_stride <- function(stride, log2_T, average) {
  if (is.null(stride)) {
    return(log2_T)
  }
  if (!isTRUE(average == "local")) {
    stop("stride only applies when averaging is local; it is incompatible ",
         "with T = 0 and T = \"global\".", call. = FALSE)
  }
  log2_stride <- log2(stride)
  if (log2_stride != round(log2_stride)) {
    stop("stride must be a power of two, got ", stride, ".", call. = FALSE)
  }
  if (log2_stride > log2_T) {
    stop("stride = ", stride, " is coarser than the averaging support T; ",
         "that would alias. Use stride <= 2^", log2_T, ".", call. = FALSE)
  }
  as.integer(log2_stride)
}

#' @noRd
check_count <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || x < 1 || x != round(x)) {
    stop(arg, " must be a single positive whole number.", call. = FALSE)
  }
  as.integer(x)
}

#' Enumerate the scattering paths of an operator
#'
#' Walks the same loops as the cascade but carries no data, so the result
#' describes the coefficients an operator will produce before any signal is
#' seen. Rows are ordered by scattering order and then by filter index, which is
#' the order [scat_transform()] returns coefficients in.
#'
#' @param op A `wavscat_op`.
#' @return A data frame with one row per path.
#' @noRd
enumerate_paths <- function(op) {
  psi1 <- op$filters$psi1
  psi2 <- op$filters$psi2

  order <- 0L
  n1 <- NA_integer_
  n2 <- NA_integer_
  j1 <- NA_integer_
  j2 <- NA_integer_

  for (i1 in seq_along(psi1)) {
    order <- c(order, 1L)
    n1 <- c(n1, i1)
    n2 <- c(n2, NA_integer_)
    j1 <- c(j1, psi1[[i1]]$j)
    j2 <- c(j2, NA_integer_)

    if (op$max_order >= 2L) {
      for (i2 in seq_along(psi2)) {
        if (psi2[[i2]]$j > psi1[[i1]]$j) {
          order <- c(order, 2L)
          n1 <- c(n1, i1)
          n2 <- c(n2, i2)
          j1 <- c(j1, psi1[[i1]]$j)
          j2 <- c(j2, psi2[[i2]]$j)
        }
      }
    }
  }

  xi1 <- ifelse(is.na(n1), NA_real_, vapply(psi1, `[[`, 0, "xi")[n1])
  xi2 <- ifelse(is.na(n2), NA_real_, vapply(psi2, `[[`, 0, "xi")[n2])
  sigma1 <- ifelse(is.na(n1), NA_real_, vapply(psi1, `[[`, 0, "sigma")[n1])
  sigma2 <- ifelse(is.na(n2), NA_real_, vapply(psi2, `[[`, 0, "sigma")[n2])

  paths <- data.frame(
    order = order, n1 = n1, n2 = n2, j1 = j1, j2 = j2,
    xi1 = xi1, xi2 = xi2, sigma1 = sigma1, sigma2 = sigma2
  )
  paths <- paths[order(paths$order, paths$n1, paths$n2, na.last = FALSE), ,
                 drop = FALSE]
  rownames(paths) <- NULL
  paths$path <- path_labels(paths)
  paths
}

#' Short human-readable label for each path
#' @noRd
path_labels <- function(paths) {
  ifelse(
    paths$order == 0L, "S0",
    ifelse(
      paths$order == 1L,
      paste0("S1_", paths$n1),
      paste0("S2_", paths$n1, "_", paths$n2)
    )
  )
}
