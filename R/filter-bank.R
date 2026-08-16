# Morlet filter bank, constructed directly in the Fourier domain.
#
# Filters are never materialised in time. Every filter is a real-valued vector
# of length N giving its Fourier transform at the frequencies of `fftfreq(N)`,
# where frequency is measured in cycles per sample so that 0.5 is Nyquist.
#
# The construction follows Anden & Mallat (2014) as implemented in Kymatio
# (kymatio/scattering1d/filter_bank.py, BSD-3-Clause). Constants and formulas
# are matched exactly so that results are numerically comparable.

#' Fourier frequencies in cycles per sample
#'
#' Equivalent to `numpy.fft.fftfreq(n)`: ascending from 0 to just below the
#' Nyquist rate, then wrapping to the negative frequencies.
#'
#' @param n Positive integer transform length.
#' @return Numeric vector of length `n`.
#' @noRd
fftfreq <- function(n) {
  n <- as.integer(n)
  half <- (n - 1L) %/% 2L
  c(0:half, if (n > 1L) seq.int(-(n %/% 2L), -1L) else integer(0)) / n
}

#' Number of periods needed to build a Morlet without wrap discontinuity
#'
#' A Gaussian of width `sigma` in frequency does not decay to zero inside a
#' single period when `sigma` is large, so evaluating it on `[0, 1)` alone
#' leaves a step at the wrap point. Choose `P` so that the Gaussian has fallen
#' below `eps` at the edges of `[1 - P, P)`, then periodise.
#'
#' @param sigma Positive bandwidth.
#' @param eps Required precision at the interval boundary.
#' @return Integer `P >= 1`; the interval `[1 - P, P)` spans `2 * P - 1` periods.
#' @noRd
adaptive_choice_P <- function(sigma, eps = 1e-7) {
  as.integer(ceiling(sqrt(-2 * sigma^2 * log(eps)) + 1))
}

#' Fourier transform of a Morlet or Gaussian filter
#'
#' In time the Morlet is `psi(t) = g_sigma(t) * (exp(1i * xi * t) - kappa)`,
#' where `g_sigma` is a Gaussian envelope and `kappa` is chosen so that `psi`
#' has exactly zero mean. Passing `xi = NULL` drops the modulation and the
#' corrective term, giving the low-pass `phi(t) = g_sigma(t)`.
#'
#' The filter is normalised so that its time-domain L1 norm is one, which makes
#' the transform non-expansive.
#'
#' @param N Transform length.
#' @param xi Centre frequency in `(0, 0.5]` cycles per sample, or `NULL` for the
#'   low-pass.
#' @param sigma Bandwidth in cycles per sample.
#' @return Real numeric vector of length `N`.
#' @noRd
morlet_1d <- function(N, xi, sigma) {
  P <- min(adaptive_choice_P(sigma), 5L)

  # Frequencies spanning 2 * P - 1 periods, laid out so that consecutive blocks
  # of length N are successive periods.
  freqs <- seq.int((1 - P) * N, P * N - 1) / N
  # With a single period, centring on [-0.5, 0.5) keeps the low-pass continuous
  # across the wrap point; the band-pass is centred on xi and does not need it.
  freqs_low <- if (P == 1L) fftfreq(N) else freqs

  low_pass_f <- fold_periods(exp(-(freqs_low^2) / (2 * sigma^2)), N)

  if (!is.null(xi) && xi != 0) {
    gabor_f <- fold_periods(exp(-(freqs - xi)^2 / (2 * sigma^2)), N)
    # Subtracting this multiple of the low-pass forces filter_f[1] == 0, i.e. a
    # wavelet with zero mean.
    kappa <- gabor_f[1L] / low_pass_f[1L]
    filter_f <- gabor_f - kappa * low_pass_f
  } else {
    filter_f <- low_pass_f
  }

  filter_f / sum(bk_modulus(bk_ifft(filter_f)))
}

#' Average a vector over successive blocks of length `N`
#'
#' Discretising in time is periodisation in frequency.
#'
#' @noRd
fold_periods <- function(v, N) {
  dim(v) <- c(N, length(v) %/% N)
  rowMeans(v)
}

#' Fourier transform of a Gaussian low-pass filter
#' @noRd
gauss_1d <- function(N, sigma) {
  morlet_1d(N, xi = NULL, sigma = sigma)
}

#' Bandwidth of a Morlet at centre frequency `xi` in a bank with `Q` per octave
#'
#' Chosen so that adjacent filters, which peak at 1, cross at height `r`. That
#' is what makes the bank tile the frequency axis without gaps.
#'
#' @param xi Centre frequency.
#' @param Q Wavelets per octave, an integer `>= 1`.
#' @param r Crossing height in `(0, 1)`.
#' @return Bandwidth `sigma`.
#' @noRd
compute_sigma_psi <- function(xi, Q, r = sqrt(0.5)) {
  factor <- 1 / 2^(1 / Q)
  term1 <- (1 - factor) / (1 + factor)
  term2 <- 1 / sqrt(2 * log(1 / r))
  xi * term1 * term2
}

#' Highest centre frequency of the bank
#'
#' Below Nyquist by enough margin that the top wavelet is not truncated.
#'
#' @noRd
compute_xi_max <- function(Q) {
  max(1 / (1 + 2^(3 / Q)), 0.35)
}

#' Largest dyadic subsampling a filter admits without aliasing
#'
#' The filter's effective upper band edge is `omega0 = |xi| + alpha * sigma`;
#' the largest `j` with `omega0 < 2^-(j + 1)` is safe.
#'
#' @param xi Centre frequency.
#' @param sigma Bandwidth.
#' @param alpha How many standard deviations to count as the band edge. Larger
#'   values mean less aliasing and less subsampling.
#' @return Non-negative integer `j`, so the filter may be subsampled by `2^j`.
#' @noRd
get_max_dyadic_subsampling <- function(xi, sigma, alpha) {
  upper_bound <- min(abs(xi) + alpha * sigma, 0.5)
  as.integer(floor(-log2(upper_bound)) - 1)
}

#' Half temporal support of a filter given in Fourier
#'
#' The smallest `N` such that truncating the filter to `[-N, N]` in time changes
#' any convolution by no more than `criterion_amplitude` relative to the input's
#' supremum. Used to decide how much to pad.
#'
#' @param h_f Real or complex vector: one filter in the Fourier domain.
#' @param criterion_amplitude Tolerated relative error.
#' @return Integer half support.
#' @noRd
compute_temporal_support <- function(h_f, criterion_amplitude = 1e-3) {
  h <- bk_ifft(h_f)
  half_support <- length(h_f) %/% 2L
  # Tail mass of |h| from each index to the middle of the support.
  l1_residual <- rev(cumsum(rev(bk_modulus(h)[seq_len(half_support)])))
  ok <- which(l1_residual <= criterion_amplitude)
  if (length(ok) == 0L) {
    warning(
      "The averaging filter is wider than the signal, so the transform cannot ",
      "be padded enough to avoid border effects. Reduce T, or use a longer ",
      "signal.",
      call. = FALSE
    )
    return(half_support)
  }
  min(ok)
}

#' Centre frequencies and bandwidths of a Morlet bank
#'
#' Two regimes, joined at an elbow frequency. Above the elbow the bank is
#' constant-Q: `xi` falls geometrically by `2^(1/Q)` per step and `sigma` follows,
#' so `xi / sigma` is fixed and the wavelets are dilations of one another. That
#' cannot continue indefinitely, because shrinking `sigma` means growing time
#' support. Below the elbow `sigma` is pinned at `sigma_min = sigma0 / 2^J` and
#' `xi` falls arithmetically instead, which bounds the time support at `2^J`
#' while still covering the frequency axis down to `2^-J`.
#'
#' @param J Log-scale of the transform; maximum wavelet support is `2^J` samples.
#' @param Q Wavelets per octave in the constant-Q region.
#' @param sigma0 Bandwidth scale; the minimum bandwidth is `sigma0 / 2^J`.
#' @param r_psi Height at which adjacent wavelets cross.
#' @return A data frame with columns `xi` and `sigma`, in descending frequency.
#' @noRd
anden_generator <- function(J, Q, sigma0, r_psi) {
  xi <- compute_xi_max(Q)
  sigma <- compute_sigma_psi(xi, Q, r = r_psi)
  sigma_min <- sigma0 / 2^J

  xis <- numeric(0)
  sigmas <- numeric(0)

  if (sigma <= sigma_min) {
    # The bank is already at its bandwidth floor: no constant-Q region at all.
    xi <- sigma
  } else {
    xis <- xi
    sigmas <- sigma
    while (sigma > sigma_min * 2^(1 / Q)) {
      xi <- xi / 2^(1 / Q)
      sigma <- sigma / 2^(1 / Q)
      xis <- c(xis, xi)
      sigmas <- c(sigmas, sigma)
    }
  }

  elbow_xi <- xi
  if (Q >= 2) {
    for (q in seq_len(Q - 1L)) {
      xi <- xi - elbow_xi / Q
      xis <- c(xis, xi)
      sigmas <- c(sigmas, sigma_min)
    }
  }

  data.frame(xi = xis, sigma = sigmas)
}

#' The bank duplicated at negative centre frequencies
#'
#' Joint time-frequency scattering needs wavelets of both signs along the
#' log-frequency axis: a filter with positive centre frequency responds to
#' spectral patterns drifting one way in time, its negative twin to the other.
#' That sign is the "spin" that lets the transform tell a rising chirp from a
#' falling one, which time scattering alone cannot do.
#'
#' @inheritParams anden_generator
#' @return A data frame of `xi` and `sigma`, positive frequencies first.
#' @noRd
spin_generator <- function(J, Q, sigma0, r_psi) {
  base <- anden_generator(J, Q, sigma0, r_psi)
  rbind(base, data.frame(xi = -base$xi, sigma = base$sigma))
}

#' Build the filter banks for a scattering transform
#'
#' Each filter is stored with its Fourier transform at several resolutions.
#' `levels[[k + 1]]` is the filter periodised for use against a signal that has
#' already been subsampled by `2^k`, which is what lets the second order operate
#' on the shortened first-order envelopes without rebuilding anything.
#'
#' @param N Padded signal length, a power of two.
#' @param J Log-scale of the transform.
#' @param Q Integer vector, one entry per layer, giving wavelets per octave.
#'   May be empty, which builds the low-pass alone.
#' @param T Support of the low-pass filter, in samples.
#' @param alpha,r_psi,sigma0 Filter bank hyperparameters.
#' @param generator Function returning the `xi` and `sigma` schedule, either
#'   [anden_generator()] or [spin_generator()].
#' @return A list with elements `phi`, one filter record, and `banks`, a list
#'   with one bank of filter records per entry of `Q`. A record has `levels`,
#'   `xi`, `sigma` and `j`.
#' @noRd
filter_factory <- function(N, J, Q, T, alpha, r_psi, sigma0,
                           generator = anden_generator) {
  log2_T <- floor(log2(T))
  max_j <- 0L
  previous_J <- 0L
  banks <- list()

  for (Q_layer in Q) {
    spec <- generator(J, Q_layer, sigma0, r_psi)
    bank <- vector("list", nrow(spec))
    for (i in seq_len(nrow(spec))) {
      xi <- spec$xi[i]
      sigma <- spec$sigma[i]
      j <- get_max_dyadic_subsampling(xi, sigma, alpha)
      # Coarser copies are only needed for filters applied after a first-order
      # subsampling, hence the cap at previous_J. The cap at j is the filter's
      # own aliasing limit; the cap at 1 + log2_T is the output stride.
      n_levels <- min(previous_J, j, 1L + log2_T)
      bank[[i]] <- list(
        levels = filter_levels(morlet_1d(N, xi, sigma), n_levels),
        xi = xi, sigma = sigma, j = j
      )
      max_j <- max(j, max_j)
    }
    previous_J <- max_j
    banks[[length(banks) + 1L]] <- bank
  }

  sigma_low <- sigma0 / T
  phi <- list(
    levels = filter_levels(gauss_1d(N, sigma_low), max(previous_J, 1L + log2_T)),
    xi = 0, sigma = sigma_low, j = log2_T, N = N
  )

  list(phi = phi, banks = banks)
}

#' Filter banks for a two-layer time scattering transform
#'
#' @inheritParams filter_factory
#' @param Q Length-2 integer vector: wavelets per octave at orders one and two.
#' @return A list with elements `phi`, `psi1` and `psi2`.
#' @noRd
scattering_filter_factory <- function(N, J, Q, T, alpha, r_psi, sigma0) {
  f <- filter_factory(N, J, Q, T, alpha, r_psi, sigma0)
  list(phi = f$phi, psi1 = f$banks[[1L]], psi2 = f$banks[[2L]])
}

#' Periodised copies of one filter at successively coarser resolutions
#'
#' @param f_full The filter at full resolution.
#' @param n_levels Number of resolutions to produce, counting the full one. A
#'   value below two yields the full-resolution filter alone.
#' @return List where element `k + 1` is the filter subsampled by `2^k`.
#' @noRd
filter_levels <- function(f_full, n_levels) {
  levels <- list(f_full)
  if (n_levels > 1L) {
    for (k in seq_len(n_levels - 1L)) {
      levels[[k + 1L]] <- bk_periodize(f_full, 2^k, reduce = "sum")
    }
  }
  levels
}

#' Fetch a filter at subsampling level `k`
#' @noRd
flevel <- function(f, k) {
  f$levels[[k + 1L]]
}
