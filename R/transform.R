# Applying an operator to signals, and the container for the result.

#' Apply a scattering operator to one or more signals
#'
#' @param op A scattering operator from [scattering_1d()].
#' @param x A numeric vector of length `op$n`, or a matrix with `op$n` rows
#'   whose columns are channels transformed independently. Column names, if
#'   present, label the channels.
#' @param ... Unused, for future extension.
#'
#' @return An object of class `wavscat_coefs`. When `op` was built with
#'   `out_type = "array"` its `coef` element is a numeric array indexed by
#'   `[path, time, channel]`; with `out_type = "list"` it is a list of matrices
#'   indexed by `[time, channel]`, one per path.
#'
#' @examples
#' op <- scattering_1d(n = 2048, J = 6, Q = 8)
#' t <- seq_len(2048) / 2048
#' x <- sin(2 * pi * 200 * t) + 0.3 * sin(2 * pi * 37 * t)
#' S <- scat_transform(op, x)
#' S
#' dim(S$coef)
#' @export
scat_transform <- function(op, x, ...) {
  if (!inherits(op, "wavscat_op")) {
    stop("op must be a scattering operator, as built by scattering_1d().",
         call. = FALSE)
  }
  x <- as_signal_matrix(x, op$n)
  channels <- colnames(x)
  n_ch <- ncol(x)

  per_channel <- lapply(seq_len(n_ch), function(ch) scat_run(op, x[, ch]))

  meta <- scat_meta(op)
  # Joint output stacked as an array keeps a frequency axis and drops order
  # zero, which has no such axis to keep.
  joint_array <- identical(op$type, "jtfs") &&
    identical(op$format, "joint") && identical(op$out_type, "array")
  if (joint_array) {
    per_channel <- lapply(per_channel, function(pc) pc[-1L])
    meta <- meta[-1L, , drop = FALSE]
  }

  coef <- if (!identical(op$out_type, "array")) {
    assemble_list(per_channel, meta$path, channels)
  } else if (joint_array) {
    assemble_joint_array(per_channel, meta$path, channels)
  } else {
    assemble_array(per_channel, meta$path, channels)
  }

  new_wavscat_coefs(coef, meta, op, channels)
}

#' Run whichever cascade the operator describes
#' @noRd
scat_run <- function(op, x) {
  if (identical(op$type, "jtfs")) {
    jtfs_unpad(op, jtfs_run(op, x))
  } else {
    scat_finalise(op, scat_cascade(op, x))
  }
}

#' Transform a signal in one step
#'
#' Builds an operator sized to `x` and applies it. Convenient for a single
#' signal; when transforming many signals of the same length, build the operator
#' once with [scattering_1d()] and reuse it, since the filter bank is the
#' expensive part.
#'
#' @param x A numeric vector, or a matrix whose columns are channels.
#' @param J,Q,T,T_sec,max_order,stride,sr,out_type Passed to [scattering_1d()].
#'   `J` defaults to the largest scale the signal supports, minus one.
#'
#' @return An object of class `wavscat_coefs`.
#'
#' @examples
#' x <- as.numeric(arima.sim(list(ar = 0.9), 2048))
#' S <- scattering(x, J = 6, Q = 8)
#' S
#' @export
scattering <- function(x, J = NULL, Q = 8, T = NULL, T_sec = NULL,
                       max_order = 2, stride = NULL, sr = NULL,
                       out_type = c("array", "list")) {
  x <- as_signal_matrix(x, n = NULL)
  n <- nrow(x)
  if (is.null(J)) {
    J <- max(floor(log2(n)) - 1L, 1L)
  }
  op <- scattering_1d(
    n = n, J = J, Q = Q, T = T, T_sec = T_sec, max_order = max_order,
    stride = stride, sr = sr, out_type = match.arg(out_type)
  )
  scat_transform(op, x)
}

#' Coerce user input to a channels-in-columns matrix
#' @noRd
as_signal_matrix <- function(x, n = NULL) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (is.null(dim(x))) {
    nm <- NULL
    x <- matrix(x, ncol = 1L, dimnames = list(NULL, nm))
  }
  if (length(dim(x)) != 2L) {
    stop("x must be a vector or a two-dimensional matrix, got an array with ",
         length(dim(x)), " dimensions.", call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop("x must be numeric.", call. = FALSE)
  }
  if (!is.null(n) && nrow(x) != n) {
    stop("This operator expects signals of length ", n, ", got ", nrow(x),
         ". Build a new operator with scattering_1d(n = ", nrow(x), ", ...).",
         call. = FALSE)
  }
  if (anyNA(x)) {
    stop("x contains missing values; scattering is not defined for them. ",
         "Interpolate or drop the affected trials first.", call. = FALSE)
  }
  if (is.null(colnames(x))) {
    colnames(x) <- if (ncol(x) == 1L) "signal" else paste0("ch", seq_len(ncol(x)))
  }
  x
}

#' @noRd
assemble_array <- function(per_channel, path_names, channels) {
  lens <- vapply(per_channel[[1L]], function(p) length(p$coef), 0L)
  if (length(unique(lens)) != 1L) {
    stop("Paths have differing lengths, so they cannot be stacked into an ",
         "array. Use out_type = \"list\".", call. = FALSE)
  }
  n_time <- lens[[1L]]
  coef <- array(
    NA_real_,
    dim = c(length(path_names), n_time, length(channels)),
    dimnames = list(path = path_names, time = NULL, channel = channels)
  )
  for (ch in seq_along(per_channel)) {
    coef[, , ch] <- matrix(
      unlist(lapply(per_channel[[ch]], `[[`, "coef"), use.names = FALSE),
      nrow = length(path_names), ncol = n_time, byrow = TRUE
    )
  }
  coef
}

#' Stack joint paths, which keep a frequency axis, into a four-way array
#' @noRd
assemble_joint_array <- function(per_channel, path_names, channels) {
  d <- dim(per_channel[[1L]][[1L]]$coef)
  shapes <- vapply(per_channel[[1L]], function(p) dim(p$coef), integer(2L))
  if (any(shapes[1L, ] != d[1L]) || any(shapes[2L, ] != d[2L])) {
    stop("Joint paths have differing shapes, so they cannot be stacked. Use ",
         "out_type = \"list\".", call. = FALSE)
  }
  coef <- array(
    NA_real_,
    dim = c(length(path_names), d[1L], d[2L], length(channels)),
    dimnames = list(path = path_names, band = NULL, time = NULL,
                    channel = channels)
  )
  for (ch in seq_along(per_channel)) {
    for (i in seq_along(path_names)) {
      coef[i, , , ch] <- per_channel[[ch]][[i]]$coef
    }
  }
  coef
}

#' @noRd
assemble_list <- function(per_channel, path_names, channels) {
  out <- lapply(seq_along(path_names), function(i) {
    coefs <- lapply(per_channel, function(pc) pc[[i]]$coef)
    first <- coefs[[1L]]
    # Channels become a trailing axis, so a plain path is a [time, channel]
    # matrix and a joint path is a [band, time, channel] array.
    d <- if (is.null(dim(first))) length(first) else dim(first)
    array(
      unlist(coefs, use.names = FALSE),
      dim = c(d, length(channels)),
      dimnames = c(rep(list(NULL), length(d)), list(channel = channels))
    )
  })
  names(out) <- path_names
  out
}

#' @noRd
new_wavscat_coefs <- function(coef, meta, op, channels) {
  structure(
    list(
      coef = coef,
      meta = meta,
      channels = channels,
      spec = list(
        type = op$type, n = op$n, J = op$J, Q = op$Q, T = op$T,
        max_order = op$max_order, average = op$average,
        log2_stride = op$log2_stride, sr = op$sr, out_type = op$out_type,
        sr_out = out_sample_rate(op)
      ),
      transforms = character(0)
    ),
    class = "wavscat_coefs"
  )
}

#' Sampling rate of the coefficients, in Hz
#' @noRd
out_sample_rate <- function(op) {
  if (is.null(op$sr)) {
    return(NA_real_)
  }
  if (identical(op$average, "global")) {
    return(NA_real_)
  }
  op$sr / 2^op$log2_stride
}
