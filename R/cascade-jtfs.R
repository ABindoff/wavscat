# The joint time-frequency scattering cascade.
#
# The shape of the computation differs from time scattering in one important
# way. Time scattering handles each first-order band independently and can
# finish one path before starting the next. Joint scattering convolves *across*
# bands, so all bands must be held together as a matrix indexed by
# [frequency band, time] before the frequential filter bank is applied. That is
# why the first order is computed breadth-first here.
#
# Passing `x = NULL` runs the cascade for its structure alone, producing the
# path table with no arithmetic. The path table and the coefficients therefore
# come from the same code and cannot drift apart.

#' Run the joint cascade over one signal, or over none
#'
#' @param op A joint scattering operator.
#' @param x Numeric vector of length `op$n`, or `NULL` to enumerate paths only.
#' @return A list of path records, sorted into output order.
#' @noRd
jtfs_run <- function(op, x = NULL) {
  dry <- is.null(x)
  phi <- op$filters$phi
  psi1 <- op$filters$psi1
  psi2 <- op$filters$psi2
  log2_stride <- op$log2_stride
  average_local <- identical(op$average, "local")

  if (!dry) {
    U_0 <- pad_reflect(x, op$pad_left, op$pad_right)
    U_0_hat <- bk_fft(U_0)
  }

  # Order zero.
  coef0 <- if (dry) {
    NULL
  } else if (average_local) {
    bk_irfft(bk_periodize(bk_cdgmm(U_0_hat, flevel(phi, 0L)), 2^log2_stride))
  } else if (identical(op$average, "global")) {
    bk_average_global(U_0)
  } else {
    U_0
  }
  zeroth <- list(list(coef = coef0, n = integer(0), order = 0L))

  # First order, breadth-first: every band is needed at once.
  U_1_hats <- vector("list", length(psi1))
  S_1_rows <- vector("list", length(psi1))
  k1s <- integer(length(psi1))
  for (n1 in seq_along(psi1)) {
    j1 <- psi1[[n1]]$j
    k1s[n1] <- if (average_local) min(j1, log2_stride) else j1
    if (dry) next
    k1 <- k1s[n1]
    U_1_hat <- bk_periodize(bk_cdgmm(U_0_hat, flevel(psi1[[n1]], 0L)), 2^k1)
    U_1_hats[[n1]] <- bk_fft(bk_modulus(bk_ifft(U_1_hat)))
    S_1_rows[[n1]] <- bk_irfft(bk_periodize(
      bk_cdgmm(U_1_hats[[n1]], flevel(phi, k1)), 2^(log2_stride - k1)
    ))
  }
  S_1 <- if (dry) NULL else do.call(rbind, S_1_rows)

  # S1 is real, so only the non-negative spins carry distinct information: the
  # negative ones would be complex conjugates.
  first <- jtfs_freq_scatter(
    list(coef = S_1, n = integer(0), n2 = NULL, j2 = NULL,
         n1_index = seq_along(psi1)),
    op, spinned = FALSE, dry = dry
  )

  second <- list()
  for (n2 in seq_along(psi2)) {
    j2 <- psi2[[n2]]$j
    rows <- list()
    n1_index <- integer(0)
    for (n1 in seq_along(psi1)) {
      if (j2 <= psi1[[n1]]$j) next
      n1_index <- c(n1_index, n1)
      if (dry) next
      k1 <- k1s[n1]
      k2 <- if (average_local) min(j2 - k1, log2_stride) else (j2 - k1)
      rows[[length(rows) + 1L]] <- bk_ifft(bk_periodize(
        bk_cdgmm(U_1_hats[[n1]], flevel(psi2[[n2]], k1)), 2^k2
      ))
    }
    if (length(n1_index) == 0L) next
    Y_2 <- if (dry) NULL else do.call(rbind, rows)
    second <- c(second, jtfs_freq_scatter(
      list(coef = Y_2, n = n2, n2 = n2, j2 = j2, n1_index = n1_index),
      op, spinned = TRUE, dry = dry
    ))
  }

  out <- c(zeroth, jtfs_finalise(c(first, second), op, dry))
  jtfs_sort(out)
}

#' Convolve along the log-frequency axis
#'
#' @param X A record holding a `[band, time]` matrix and the bands it covers.
#' @param spinned Whether the input is complex. Real input, which is the first
#'   order, only needs the non-negative spins; the negative ones are conjugates
#'   and carry nothing new.
#' @return A list of path records, one per frequential filter.
#' @noRd
jtfs_freq_scatter <- function(X, op, spinned, dry) {
  psis_fr <- op$filters_fr$psis
  n_fr_pad <- op$filters_fr$phi$N
  log2_stride_fr <- op$log2_stride_fr
  average_local_fr <- identical(op$average_fr, "local")
  n1_max <- length(X$n1_index)

  if (!dry) {
    X_hat <- bk_fft_freq(rbind(
      X$coef,
      matrix(0, nrow = n_fr_pad - n1_max, ncol = ncol(X$coef))
    ))
  }

  out <- list()
  for (n_fr in seq_along(psis_fr)) {
    psi <- psis_fr[[n_fr]]
    if (!spinned && psi$xi < 0) next
    j_fr <- psi$j
    k_fr <- if (average_local_fr) min(j_fr, log2_stride_fr) else j_fr
    coef <- if (dry) {
      NULL
    } else {
      bk_ifft_freq(bk_periodize_axis(
        bk_cdgmm(X_hat, psi$levels[[1L]]), 2^k_fr, axis = 1L
      ))
    }
    out[[length(out) + 1L]] <- list(
      coef = coef, n = c(X$n, n_fr), n2 = X$n2, j2 = X$j2,
      n_fr = n_fr, j_fr = j_fr, spin = sign(psi$xi),
      n1_index = X$n1_index, n1_max = n1_max,
      n1_stride = as.integer(2^j_fr)
    )
  }
  out
}

#' Take moduli, average in time and frequency, and lay the paths out
#' @noRd
jtfs_finalise <- function(paths, op, dry) {
  phi <- op$filters$phi
  phi_fr <- op$filters_fr$phi
  average <- op$average
  average_fr <- op$average_fr
  keep_freq_axis <- identical(op$format, "joint")

  out <- list()
  for (p in paths) {
    if (!dry) {
      p$coef <- bk_modulus(p$coef)
    }

    # Temporal averaging. First-order paths descend from S1, which was already
    # averaged, so only the second order needs it here.
    if (identical(average, "global")) {
      if (!dry) p$coef <- matrix(rowSums(p$coef), ncol = 1L)
    } else if (identical(average, "local") && length(p$n) > 1L) {
      if (!dry) p$coef <- jtfs_time_average(p, phi, op$log2_stride)
    }

    # Frequential averaging. The unspun path came from a frequential low-pass
    # already, so averaging it again would be redundant.
    if (!isFALSE(average_fr) && p$spin != 0) {
      if (identical(average_fr, "global")) {
        if (!dry) p$coef <- matrix(colSums(p$coef), nrow = 1L)
        p$n1_stride <- as.integer(p$n1_max)
      } else {
        if (!dry) p$coef <- jtfs_freq_average(p, phi_fr, op$log2_F,
                                              op$log2_stride_fr)
        p$n1_stride <- as.integer(2^op$log2_stride_fr)
      }
    } else if (isFALSE(average_fr)) {
      p$n1_stride <- as.integer(2^max(p$j_fr, 0L))
    }

    if (!keep_freq_axis || !identical(op$out_type, "array")) {
      n_keep <- 1L + p$n1_max %/% p$n1_stride
      if (!dry) p$coef <- p$coef[seq_len(min(n_keep, nrow(p$coef))), , drop = FALSE]
    }

    if (keep_freq_axis) {
      p$order <- length(p$n)
      out[[length(out) + 1L]] <- p
    } else {
      out <- c(out, jtfs_split_bands(p, dry))
    }
  }
  out
}

#' Average a joint path along time
#'
#' Done in a transposed, time-major layout. `mvfft` transforms columns, so
#' putting time down the columns lets the whole chain run without further
#' transposes and lets the filter recycle naturally; only the two transposes at
#' the ends remain.
#'
#' @noRd
jtfs_time_average <- function(p, phi, log2_stride) {
  k_in <- p$j2
  k_J <- log2_stride - k_in
  if (k_J < 0L) {
    stop("The averaging support T is shorter than the widest second-order ",
         "wavelet, which joint scattering cannot accommodate. Increase T, or ",
         "reduce J.", call. = FALSE)
  }
  Mt <- t(p$coef)
  H <- stats::mvfft(Mt) * flevel(phi, k_in)
  H <- bk_periodize_axis(H, 2^k_J, axis = 1L)
  t(Re(stats::mvfft(H, inverse = TRUE) / nrow(H)))
}

#' @noRd
jtfs_freq_average <- function(p, phi_fr, log2_F, log2_stride_fr) {
  k_in <- min(p$j_fr, log2_F)
  k_J <- log2_stride_fr - k_in
  Re(bk_ifft_freq(bk_periodize_axis(
    bk_cdgmm(bk_fft_freq(p$coef), flevel(phi_fr, k_in)), 2^k_J, axis = 1L
  )))
}

#' Split a joint path into one series per frequency band
#'
#' Used by `format = "time"`, which trades the 2D structure for output shaped
#' like time scattering.
#'
#' @noRd
jtfs_split_bands <- function(p, dry) {
  starts <- as.integer(seq.int(0L, p$n1_max - 1L, by = p$n1_stride))
  lapply(seq_along(starts), function(i) {
    list(
      coef = if (dry) NULL else p$coef[i, ],
      n = c(starts[i], p$n),
      # n1 is reported one-based, and refers to the first band of the group
      # this coefficient summarises.
      n1 = p$n1_index[starts[i] + 1L],
      n2 = p$n2, n_fr = p$n_fr, j2 = p$j2, j_fr = p$j_fr, spin = p$spin,
      # The split adds n1 to the front of the path, so the post-split path
      # length is one more than p$n and the order is length(p$n).
      order = length(p$n)
    )
  })
}

#' Trim the joint cascade output back to the extent of the input signal
#'
#' @param op A joint scattering operator.
#' @param paths Output of [jtfs_run()].
#' @return The same list, with `coef` unpadded along time.
#' @noRd
jtfs_unpad <- function(op, paths) {
  average <- op$average
  if (identical(average, "global")) {
    return(paths)
  }
  for (i in seq_along(paths)) {
    p <- paths[[i]]
    res <- if (p$order == 0L) {
      if (isFALSE(average)) 0L else op$log2_stride
    } else if (isFALSE(average) && length(p$n) > 1L) {
      # Unaveraged second-order paths sit at the resolution of their widest
      # wavelet; first-order ones were already brought to the output stride.
      max(p$j2 %||% -1L, 0L)
    } else {
      max(op$log2_stride, 0L)
    }
    paths[[i]]$coef <- if (is.matrix(p$coef)) {
      p$coef[, unpad_index(op$borders, res), drop = FALSE]
    } else {
      unpad(p$coef, op$borders, res)
    }
  }
  paths
}

#' Sort paths into output order: by order, then by filter indices
#' @noRd
jtfs_sort <- function(paths) {
  lens <- vapply(paths, function(p) length(p$n), 0L)
  max_len <- max(lens)
  keys <- lapply(seq_len(max_len), function(k) {
    vapply(paths, function(p) if (length(p$n) >= k) as.numeric(p$n[k]) else -1, 0)
  })
  paths[do.call(order, c(list(lens), keys))]
}
