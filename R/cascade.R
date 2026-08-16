# The scattering cascade.
#
# Wavelet convolution, complex modulus, low-pass averaging, iterated. Everything
# stays in the Fourier domain so that subsampling is free: periodising a
# spectrum is exactly subsampling the signal it represents, which keeps every
# inverse transform as short as the bandwidth allows. That is what makes the
# transform affordable without compiled code.

#' Run the cascade on one padded signal
#'
#' @param op A `wavscat_op`.
#' @param x Numeric vector of length `op$n`.
#' @return A list of path records with elements `coef`, `order`, `n1`, `n2`,
#'   `j1`, `j2`, in generation order.
#' @noRd
scat_cascade <- function(op, x) {
  phi <- op$filters$phi
  psi1 <- op$filters$psi1
  psi2 <- op$filters$psi2
  log2_stride <- op$log2_stride
  average_local <- identical(op$average, "local")
  second_order <- op$max_order >= 2L

  out <- vector("list", nrow(op$paths))
  k <- 0L
  emit <- function(coef, order, n1 = NA_integer_, n2 = NA_integer_,
                   j1 = NA_integer_, j2 = NA_integer_) {
    k <<- k + 1L
    out[[k]] <<- list(coef = coef, order = order, n1 = n1, n2 = n2,
                      j1 = j1, j2 = j2)
  }

  U_0 <- pad_reflect(x, op$pad_left, op$pad_right)
  U_0_hat <- bk_fft(U_0)

  # Order zero: the local mean of the signal.
  if (average_local) {
    S_0 <- bk_irfft(bk_periodize(bk_cdgmm(U_0_hat, flevel(phi, 0L)), 2^log2_stride))
    emit(S_0, 0L)
  } else {
    emit(U_0, 0L)
  }

  for (n1 in seq_along(psi1)) {
    j1 <- psi1[[n1]]$j
    # Subsample as far as the wavelet's bandwidth allows, but no further than
    # the output stride, since anything finer would be averaged away anyway.
    k1 <- if (average_local) min(j1, log2_stride) else j1

    U_1_hat <- bk_periodize(bk_cdgmm(U_0_hat, flevel(psi1[[n1]], 0L)), 2^k1)
    U_1_m <- bk_modulus(bk_ifft(U_1_hat))

    # The modulus is a smooth envelope, so its spectrum is reused both for the
    # averaging that produces S1 and for the second-order convolutions.
    if (average_local || second_order) {
      U_1_hat <- bk_fft(U_1_m)
    }

    if (average_local) {
      S_1 <- bk_irfft(bk_periodize(
        bk_cdgmm(U_1_hat, flevel(phi, k1)), 2^max(log2_stride - k1, 0L)
      ))
      emit(S_1, 1L, n1 = n1, j1 = j1)
    } else {
      emit(U_1_m, 1L, n1 = n1, j1 = j1)
    }

    if (!second_order) next

    for (n2 in seq_along(psi2)) {
      j2 <- psi2[[n2]]$j
      # |x * psi_j1| is an envelope band-limited to roughly the bandwidth of
      # psi_j1, so second-order wavelets at or above that frequency capture
      # nothing. Skipping them is what keeps the path count near J^2 Q / 2
      # rather than (J Q)^2.
      if (j2 <= j1) next

      sub2_adj <- if (average_local) min(j2, log2_stride) else j2
      k2 <- max(sub2_adj - k1, 0L)

      U_2_hat <- bk_periodize(bk_cdgmm(U_1_hat, flevel(psi2[[n2]], k1)), 2^k2)
      U_2_m <- bk_modulus(bk_ifft(U_2_hat))

      if (average_local) {
        S_2 <- bk_irfft(bk_periodize(
          bk_cdgmm(bk_fft(U_2_m), flevel(phi, k1 + k2)),
          2^max(log2_stride - sub2_adj, 0L)
        ))
        emit(S_2, 2L, n1 = n1, n2 = n2, j1 = j1, j2 = j2)
      } else {
        emit(U_2_m, 2L, n1 = n1, n2 = n2, j1 = j1, j2 = j2)
      }
    }
  }

  out[seq_len(k)]
}

#' Trim the cascade output back to the extent of the input signal
#'
#' @param op A `wavscat_op`.
#' @param paths Output of [scat_cascade()].
#' @return The same list, with `coef` unpadded or globally averaged, and
#'   reordered to match `op$paths`.
#' @noRd
scat_finalise <- function(op, paths) {
  average_local <- identical(op$average, "local")
  average_global <- identical(op$average, "global")

  for (i in seq_along(paths)) {
    p <- paths[[i]]
    if (average_global) {
      # Note that this sums over the padded signal, matching Kymatio.
      paths[[i]]$coef <- bk_average_global(p$coef)
    } else {
      res <- if (average_local) {
        op$log2_stride
      } else if (p$order > 0L) {
        if (p$order == 2L) p$j2 else p$j1
      } else {
        0L
      }
      paths[[i]]$coef <- unpad(p$coef, op$borders, res)
    }
  }

  key <- order(
    vapply(paths, `[[`, 0L, "order"),
    vapply(paths, function(p) if (is.na(p$n1)) -1L else p$n1, 0L),
    vapply(paths, function(p) if (is.na(p$n2)) -1L else p$n2, 0L)
  )
  paths[key]
}
