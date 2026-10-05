# Turning coefficients into model-ready features.

#' Log-compress scattering coefficients
#'
#' Coefficients of order one and above are magnitudes, so their distribution
#' across paths spans several orders of magnitude and is strongly right-skewed.
#' Compressing with `log(1 + x / eps)` brings them onto a comparable scale,
#' which is what most classifiers want, and turns the multiplicative effect of
#' overall signal amplitude into an additive offset.
#'
#' The order-zero coefficient is a local mean rather than a magnitude and may be
#' negative, so the compression applied is `sign(x) * log(1 + |x| / eps)`. This
#' agrees exactly with `log(1 + x / eps)` wherever `x` is non-negative, and
#' stays monotone and finite where it is not.
#'
#' @section Choosing eps:
#'
#' `eps` sets the value below which coefficients are treated as noise: well
#' above `eps` the transform is essentially `log(x / eps)`, and well below it is
#' essentially linear. It must therefore sit somewhere near the noise floor of
#' your coefficients, which depends on how your signals are scaled.
#'
#' **Use one `eps` for a whole dataset, and derive it from training data only.**
#' Choosing it separately per recording, or from the full dataset, leaks
#' information between train and test. [scat_eps()] computes a reasonable value
#' from a set of training coefficients; [step_scattering()] does this
#' automatically at preparation time.
#'
#' @param x A `wavscat_coefs` object.
#' @param eps Positive scalar; the scale at which compression begins.
#'
#' @return A `wavscat_coefs` object with compressed coefficients.
#'
#' @seealso [scat_eps()], [scat_renorm()]
#'
#' @examples
#' S <- scattering(as.numeric(arima.sim(list(ar = 0.9), 2048)), J = 6, Q = 8)
#' scat_log(S, eps = scat_eps(S))
#' @export
scat_log <- function(x, eps = 1e-6) {
  check_coefs(x)
  if (!is.numeric(eps) || length(eps) != 1L || eps <= 0) {
    stop("eps must be a single positive number.", call. = FALSE)
  }
  if ("log" %in% x$transforms) {
    warning("These coefficients have already been log-compressed.", call. = FALSE)
  }
  x$coef <- map_coef(x$coef, function(v) sign(v) * log1p(abs(v) / eps))
  x$transforms <- c(x$transforms, "log")
  x
}

#' A default compression scale for a set of coefficients
#'
#' Returns a small quantile of the non-zero first-order coefficients, which puts
#' the elbow of [scat_log()] near the bottom of the observed range rather than
#' at an arbitrary absolute value. Compute this once on training data and pass
#' the result to `scat_log()` for every split.
#'
#' @param x A `wavscat_coefs` object, or a list of them.
#' @param quantile Quantile of the first-order coefficients to use.
#' @return A single positive number.
#' @examples
#' S <- scattering(as.numeric(arima.sim(list(ar = 0.9), 2048)), J = 6, Q = 8)
#' scat_eps(S)
#' @export
scat_eps <- function(x, quantile = 0.01) {
  if (inherits(x, "wavscat_coefs")) {
    x <- list(x)
  }
  vals <- unlist(lapply(x, function(s) {
    check_coefs(s)
    keep <- s$meta$order >= 1L
    if (is.array(s$coef)) as.numeric(s$coef[keep, , , drop = FALSE])
    else unlist(s$coef[keep], use.names = FALSE)
  }), use.names = FALSE)
  vals <- vals[is.finite(vals) & vals > 0]
  if (length(vals) == 0L) {
    return(1e-6)
  }
  max(stats::quantile(vals, quantile, names = FALSE), .Machine$double.eps)
}

#' Renormalise second-order coefficients by their first-order parent
#'
#' Divides every order-2 coefficient by the order-1 coefficient it descends
#' from, path `n1`. Without this, an order-2 coefficient mostly reports how much
#' energy the band `n1` carried; after it, the coefficient reports the *shape*
#' of that band's amplitude modulation, largely independent of how loud the band
#' was. Anden and Mallat (2014) report a substantial reduction in classification
#' error from this step.
#'
#' Apply it before [scat_log()], on raw magnitudes.
#'
#' Time scattering only, for now. Joint time-frequency coefficients are refused,
#' because a joint first-order path is filtered along frequency and so is not
#' the parent of any second-order path.
#'
#' @param x A `wavscat_coefs` object.
#' @param eps Small constant added to the denominator to keep silent bands from
#'   producing enormous ratios. Note the trade-off: with `eps = 0` the
#'   renormalised coefficients are exactly invariant to the amplitude of the
#'   recording, but a band carrying no energy divides by zero. Any positive
#'   `eps` restores stability and breaks that invariance, mildly for loud bands
#'   and completely for silent ones, which is the intended behaviour.
#'
#' @return A `wavscat_coefs` object.
#'
#' @examples
#' S <- scattering(as.numeric(arima.sim(list(ar = 0.9), 2048)), J = 6, Q = 8)
#' S <- scat_renorm(S)
#' S <- scat_log(S, eps = scat_eps(S))
#' @export
scat_renorm <- function(x, eps = 1e-12) {
  check_coefs(x)
  if (identical(x$spec$type, "jtfs")) {
    # A joint first-order path is itself filtered along frequency, so there
    # are several per band and none is the parent of a second-order path. The
    # right denominator, S1 of the bands a path spans through the same
    # frequential low-pass, is defined in wavscat-core and arrives with the
    # Rust engine; until then, refuse rather than divide by the wrong thing.
    stop("scat_renorm() does not yet support joint time-frequency ",
         "coefficients: their first-order paths are filtered along frequency, ",
         "so none is the parent of a second-order path. Use scat_log() alone ",
         "for now.", call. = FALSE)
  }
  if ("log" %in% x$transforms) {
    stop("Renormalise before log-compressing: dividing log magnitudes does not ",
         "give the intended ratio.", call. = FALSE)
  }
  if ("renorm" %in% x$transforms) {
    warning("These coefficients have already been renormalised.", call. = FALSE)
  }
  meta <- x$meta
  o2 <- which(meta$order == 2L)
  if (length(o2) == 0L) {
    x$transforms <- c(x$transforms, "renorm")
    return(x)
  }
  # Each order-2 path has exactly one order-1 parent, sharing its n1.
  parent_of <- match(meta$n1[o2], meta$n1[meta$order == 1L])
  parent_row <- which(meta$order == 1L)[parent_of]

  if (is.array(x$coef)) {
    x$coef[o2, , ] <- x$coef[o2, , , drop = FALSE] /
      (x$coef[parent_row, , , drop = FALSE] + eps)
  } else {
    for (k in seq_along(o2)) {
      x$coef[[o2[k]]] <- x$coef[[o2[k]]] / (x$coef[[parent_row[k]]] + eps)
    }
  }
  x$transforms <- c(x$transforms, "renorm")
  x
}

#' Extract a tidy feature table
#'
#' One row per channel, one column per feature, with names that decode back to
#' the path they came from. Feed the result straight into a modelling workflow.
#'
#' @param x A `wavscat_coefs` object with array storage.
#' @param summary How to reduce the time axis of each path. `"mean"`, the
#'   default, gives one feature per path. `"none"` keeps every time point, so
#'   the number of features is paths times time points. `"max"`, `"sd"` and
#'   `"median"` are also available; `"sd"` measures how much a path fluctuates
#'   over the recording, which for a tapping trace is a measure of rhythm
#'   variability.
#' @param orders Which scattering orders to keep. Defaults to all of them.
#' @param prefix Optional string prepended to every feature name, useful when
#'   binding features from several operators side by side.
#'
#' @return A tibble whose first column, `channel`, names the channel, followed
#'   by the feature columns.
#'
#' @examples
#' set.seed(1)
#' x <- cbind(thumb = rnorm(2048), index = rnorm(2048))
#' S <- scattering(x, J = 6, Q = 8)
#' f <- scat_features(scat_log(S, eps = scat_eps(S)))
#' dim(f)
#' f[, 1:5]
#'
#' # Keep the time axis, first order only.
#' scat_features(S, summary = "none", orders = 1)
#' @export
scat_features <- function(x, summary = c("mean", "none", "max", "sd", "median"),
                          orders = NULL, prefix = NULL) {
  check_coefs(x)
  summary <- match.arg(summary)
  if (!is.array(x$coef)) {
    stop("scat_features() needs array storage; rebuild the operator with ",
         "out_type = \"array\".", call. = FALSE)
  }
  if (length(dim(x$coef)) == 4L) {
    stop("Joint coefficients keep a frequency axis, which does not flatten to ",
         "one feature per path. Build the operator with format = \"time\".",
         call. = FALSE)
  }

  keep <- if (is.null(orders)) rep(TRUE, nrow(x$meta)) else x$meta$order %in% orders
  if (!any(keep)) {
    stop("No paths of order ", paste(orders, collapse = ", "), " to extract.",
         call. = FALSE)
  }
  coef <- x$coef[keep, , , drop = FALSE]
  labels <- x$meta$path[keep]

  if (summary == "none") {
    n_time <- dim(coef)[2L]
    # [path, time, channel] to [channel, path * time], with time varying fastest
    # so that each path's time course stays contiguous.
    values <- t(matrix(aperm(coef, c(2L, 1L, 3L)),
                       nrow = length(labels) * n_time,
                       ncol = dim(coef)[3L]))
    names_out <- paste0(rep(labels, each = n_time), "_t", seq_len(n_time))
  } else {
    fn <- switch(summary,
      mean = mean, max = max, sd = stats::sd, median = stats::median
    )
    values <- t(apply(coef, c(1L, 3L), fn))
    if (is.null(dim(values))) {
      values <- matrix(values, nrow = dim(coef)[3L])
    }
    names_out <- labels
  }

  if (!is.null(prefix)) {
    names_out <- paste0(prefix, names_out)
  }
  colnames(values) <- names_out

  tibble::as_tibble(cbind(
    data.frame(channel = x$channels, stringsAsFactors = FALSE),
    as.data.frame(values)
  ))
}

#' @noRd
check_coefs <- function(x) {
  if (!inherits(x, "wavscat_coefs")) {
    stop("x must be a wavscat_coefs object, as returned by scat_transform().",
         call. = FALSE)
  }
  invisible(x)
}

#' Apply a function to coefficients held either as an array or as a list
#' @noRd
map_coef <- function(coef, fn) {
  if (is.array(coef)) {
    out <- fn(coef)
    dim(out) <- dim(coef)
    dimnames(out) <- dimnames(coef)
    out
  } else {
    lapply(coef, fn)
  }
}
