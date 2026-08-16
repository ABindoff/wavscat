# Construction of a joint time-frequency scattering operator.

#' Create a joint time-frequency scattering operator
#'
#' Time scattering treats each frequency band separately: it records how the
#' energy in each band fluctuates, but not how bands move together. Joint
#' time-frequency scattering adds a second filter bank running along the
#' log-frequency axis, so the second order convolves jointly in time and in
#' frequency.
#'
#' What that buys is sensitivity to *spectrotemporal* structure. The frequential
#' wavelets come in two spins, positive and negative, which respond to spectral
#' patterns drifting up and down in frequency respectively. A rising formant
#' transition and a falling one have identical time-scattering coefficients and
#' opposite spins here. For speech and DDK recordings, where the informative
#' cues are exactly these formant movements and the coordination of energy
#' across bands, that difference matters.
#'
#' The cost is more coefficients and more computation, so prefer
#' [scattering_1d()] when the signal has no meaningful structure across
#' frequency, which is usually the case for a single kinematic trace.
#'
#' @param n Length of the signals to be transformed, in samples.
#' @param J Log-scale of the temporal transform, as in [scattering_1d()].
#' @param J_fr Log-scale of the frequential transform. The widest frequential
#'   wavelet spans `2^J_fr` first-order bands, so with `Q = 8` a `J_fr` of 3
#'   reaches across one octave.
#' @param Q Wavelets per octave in time; one number, or `c(Q1, Q2)`.
#' @param Q_fr Wavelets per octave along log-frequency. One is usually right;
#'   the frequency axis is short and rarely supports finer resolution.
#' @param T Temporal averaging support in samples, `NULL` for `2^J`, or
#'   `"global"`.
#'
#'   Avoid `"global"` here in particular. It sums over the reflection-padded
#'   signal, and the mirror image of a rising sweep is a falling one, so the two
#'   spins receive equal energy and the spin asymmetry that motivates this
#'   transform is cancelled almost exactly. Keep `T` local and reduce the time
#'   axis afterwards with `scat_features(summary = "mean")`.
#' @param T_sec Temporal averaging support in seconds; requires `sr`.
#' @param F Frequential averaging support, measured in first-order bands.
#'   `NULL` uses `2^J_fr`, `"global"` collapses the frequency axis entirely, and
#'   `0` leaves it unaveraged.
#' @param stride,stride_fr Output subsampling factors in time and in frequency,
#'   each a power of two.
#' @param format `"time"` returns one coefficient series per
#'   `(band, modulation, spin)` path, which stacks into the same shape as
#'   [scattering_1d()] output and feeds the same downstream functions.
#'   `"joint"` keeps the frequency axis intact, giving a 2D image per path,
#'   which is what a convolutional model would want.
#' @param sr Sampling rate in Hz, used only for reporting.
#' @param out_type `"array"` or `"list"`.
#'
#' @return An object of class `wavscat_op`.
#'
#' @references
#' Anden, J., Lostanlen, V. and Mallat, S. (2019). Joint time-frequency
#' scattering. *IEEE Transactions on Signal Processing*, 67(14), 3704-3718.
#' \doi{10.1109/TSP.2019.2918992}
#'
#' @seealso [scattering_1d()] for time scattering.
#'
#' @examples
#' op <- scattering_jtfs(n = 4096, J = 6, J_fr = 3, Q = 8, sr = 16000)
#' op
#'
#' # Spins are balanced: as many rising as falling, plus the unspun paths.
#' table(scat_meta(op)$spin)
#' @export
scattering_jtfs <- function(n, J, J_fr = 3, Q = 8, Q_fr = 1, T = NULL,
                            T_sec = NULL, F = NULL, stride = NULL,
                            stride_fr = NULL, format = c("time", "joint"),
                            sr = NULL, out_type = c("array", "list")) {
  format <- match.arg(format)
  out_type <- match.arg(out_type)

  # The temporal half of the operator is exactly a time scattering operator, so
  # build one and extend it rather than duplicating the validation.
  op <- scattering_1d(
    n = n, J = J, Q = Q, T = T, T_sec = T_sec, max_order = 2L,
    stride = stride, sr = sr, out_type = out_type
  )

  J_fr <- check_count(J_fr, "J_fr")
  Q_fr <- check_count(Q_fr, "Q_fr")

  if (format == "joint" && out_type == "array" && identical(F, 0)) {
    stop("format = \"joint\" with out_type = \"array\" needs frequential ",
         "averaging, because unaveraged paths have differing numbers of bands. ",
         "Set F to a positive value, or use out_type = \"list\".", call. = FALSE)
  }

  # Kymatio sizes the frequency axis from the nominal band count rather than the
  # realised one, so that F and the padding do not shift when the bank changes
  # by a filter at the edge.
  n_input_fr <- (op$J + 1L) * op$Q[1L]
  if (2^J_fr > n_input_fr) {
    stop(
      "J_fr = ", J_fr, " asks for a frequential wavelet spanning 2^", J_fr,
      " = ", 2^J_fr, " bands, but this operator has only ", n_input_fr,
      " first-order bands. Use J_fr <= ", floor(log2(n_input_fr)),
      ", or raise J or Q to widen the frequency axis.",
      call. = FALSE
    )
  }
  parsed_fr <- parse_T(F, J_fr, n_input_fr, arg = "F")
  F <- parsed_fr$T
  average_fr <- parsed_fr$average
  log2_F <- as.integer(floor(log2(F)))
  log2_stride_fr <- parse_stride_fr(stride_fr, log2_F, average_fr)

  # Enough room that the frequential convolution does not wrap, and a multiple
  # of every subsampling factor it will use.
  min_to_pad_fr <- 8 * min(F, 2^J_fr)
  K_fr <- max(J_fr, 0L)
  n_padded_fr <- as.integer(
    ((n_input_fr + min_to_pad_fr) %/% 2^K_fr) * 2^K_fr
  )

  phi_fr <- filter_factory(
    n_padded_fr, J_fr, integer(0), F, op$alpha, op$r_psi, op$sigma0,
    generator = spin_generator
  )$phi
  spun <- filter_factory(
    n_padded_fr, J_fr, Q_fr, 2^J_fr, op$alpha, op$r_psi, op$sigma0,
    generator = spin_generator
  )
  # Index 1 is the unspun low-pass, which carries the paths that average across
  # frequency rather than differentiating along it.
  psis_fr <- c(list(spun$phi), spun$banks[[1L]])

  aliased <- vapply(spun$banks[[1L]],
                    function(p) abs(p$xi) >= 0.5 / 2^p$j, TRUE)
  if (any(aliased)) {
    stop("The frequential filter bank would alias with J_fr = ", J_fr,
         " and Q_fr = ", Q_fr, ". Reduce J_fr.", call. = FALSE)
  }

  op$type <- "jtfs"
  op$J_fr <- J_fr
  op$Q_fr <- Q_fr
  op$F <- F
  op$average_fr <- average_fr
  op$log2_F <- log2_F
  op$log2_stride_fr <- log2_stride_fr
  op$n_input_fr <- n_input_fr
  op$n_padded_fr <- n_padded_fr
  op$format <- format
  op$filters_fr <- list(phi = phi_fr, psis = psis_fr)

  op$paths <- jtfs_enumerate_paths(op)
  op
}

#' @noRd
parse_stride_fr <- function(stride_fr, log2_F, average_fr) {
  if (is.null(stride_fr)) {
    return(log2_F)
  }
  if (!isTRUE(average_fr == "local")) {
    stop("stride_fr only applies when frequential averaging is local; it is ",
         "incompatible with F = 0 and F = \"global\".", call. = FALSE)
  }
  v <- log2(stride_fr)
  if (v != round(v)) {
    stop("stride_fr must be a power of two, got ", stride_fr, ".", call. = FALSE)
  }
  if (v > log2_F) {
    stop("stride_fr = ", stride_fr, " is coarser than the frequential ",
         "averaging support F. Use stride_fr <= 2^", log2_F, ".", call. = FALSE)
  }
  as.integer(v)
}

#' Enumerate the paths of a joint operator
#'
#' Runs the cascade with no signal attached, so the path table is produced by
#' exactly the code that will later produce the coefficients. Deriving it
#' independently would risk the two drifting apart.
#'
#' @noRd
jtfs_enumerate_paths <- function(op) {
  paths <- jtfs_run(op, x = NULL)
  psi1 <- op$filters$psi1
  psi2 <- op$filters$psi2
  psis_fr <- op$filters_fr$psis

  pick <- function(bank, idx, field) {
    ifelse(is.na(idx), NA_real_,
           vapply(bank, function(p) as.numeric(p[[field]]), 0)[idx])
  }

  n1 <- vapply(paths, function(p) as.integer(p$n1 %||% NA_integer_), 0L)
  n2 <- vapply(paths, function(p) as.integer(p$n2 %||% NA_integer_), 0L)
  n_fr <- vapply(paths, function(p) as.integer(p$n_fr %||% NA_integer_), 0L)

  out <- data.frame(
    order = vapply(paths, `[[`, 0L, "order"),
    n1 = n1, n2 = n2, n_fr = n_fr,
    j1 = pick(psi1, n1, "j"), j2 = pick(psi2, n2, "j"),
    j_fr = pick(psis_fr, n_fr, "j"),
    xi1 = pick(psi1, n1, "xi"), xi2 = pick(psi2, n2, "xi"),
    xi_fr = pick(psis_fr, n_fr, "xi"),
    sigma1 = pick(psi1, n1, "sigma"), sigma2 = pick(psi2, n2, "sigma"),
    sigma_fr = pick(psis_fr, n_fr, "sigma"),
    spin = vapply(paths, function(p) as.numeric(p$spin %||% NA_real_), 0)
  )
  out$path <- jtfs_path_labels(out, op$format)
  out
}

#' @noRd
jtfs_path_labels <- function(paths, format) {
  # u for a wavelet spinning up in frequency, d for down, n for the unspun
  # low-pass path.
  spin_tag <- ifelse(
    is.na(paths$spin), "",
    ifelse(paths$spin > 0, "u", ifelse(paths$spin < 0, "d", "n"))
  )
  if (format == "joint") {
    ifelse(
      paths$order == 0L, "S0",
      ifelse(paths$order == 1L,
             paste0("J1_", paths$n_fr, spin_tag),
             paste0("J2_", paths$n2, "_", paths$n_fr, spin_tag))
    )
  } else {
    ifelse(
      paths$order == 0L, "S0",
      ifelse(paths$order == 1L,
             paste0("J1_", paths$n1, "_", paths$n_fr, spin_tag),
             paste0("J2_", paths$n1, "_", paths$n2, "_", paths$n_fr, spin_tag))
    )
  }
}

#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x
