# Presets and signal preparation.
#
# The core is signal-agnostic. What lives here is the small amount of knowledge
# that is specific to a kind of recording: sensible scales for webcam-rate
# kinematics versus speech audio, and the two preparation steps that behavioural
# recordings almost always need before a fixed-length operator can be applied.

PRESETS <- list(
  kinematic = list(
    max_scale_sec = 4, T_sec = 1, Q = c(8L, 1L),
    note = "webcam-rate movement, such as MediaPipe landmark traces at ~30 fps"
  ),
  speech = list(
    max_scale_sec = 0.25, T_sec = 0.025, Q = c(8L, 2L),
    note = "speech and diadochokinetic audio at 16 kHz or above"
  ),
  audio = list(
    max_scale_sec = 0.19, T_sec = 0.025, Q = c(8L, 2L),
    note = "general audio at 44.1 kHz"
  )
)

#' Scattering operators with defaults for a kind of recording
#'
#' Wraps [scattering_1d()] with parameter choices that suit a class of signal.
#' The choices are ordinary arguments, so override any of them by naming it.
#'
#' @section What the presets choose, and why:
#'
#' \describe{
#'   \item{`"kinematic"`}{For movement traces sampled at video rates. The largest
#'     wavelet spans about 4 seconds, which reaches well below a finger-tapping
#'     rate of 2 to 7 Hz and so captures drift and fatigue across a trial as well
#'     as the taps themselves. `Q = c(8, 1)` resolves tapping rate finely at the
#'     first order while keeping the second order broad enough to localise
#'     individual taps. Averaging over 1 second keeps within-trial structure:
#'     the default `T` of `2^J` would average over the full 4 seconds and flatten
#'     exactly the variation of interest.}
#'   \item{`"speech"`}{For speech and DDK audio. The largest wavelet spans 250 ms,
#'     so second-order paths reach modulation rates down to about 4 Hz, which
#'     covers the syllable rate. `Q = c(8, 2)` follows Anden and Mallat, whose
#'     `Q1 = 8` gives the frequency resolution needed to separate nearby
#'     harmonics and whose `Q2 = 2` resolves short modulated structures. The
#'     25 ms averaging retains the syllable train rather than smoothing it away.}
#'   \item{`"audio"`}{As `"speech"`, with a slightly shorter maximum scale suited
#'     to 44.1 kHz material.}
#' }
#'
#' These are starting points, not recommendations for your data. Plot the
#' operator with [plot.wavscat_op()] and check [scat_energy()] on a few real
#' recordings before settling.
#'
#' @param kind One of `"kinematic"`, `"speech"` or `"audio"`.
#' @param n Signal length in samples.
#' @param sr Sampling rate in Hz.
#' @param ... Any argument of [scattering_1d()], overriding the preset.
#'
#' @return A `wavscat_op`.
#'
#' @examples
#' # A 60 second tapping trial at 30 fps.
#' op <- scat_preset("kinematic", n = 1800, sr = 30)
#' op
#'
#' # Two seconds of DDK audio at 16 kHz.
#' scat_preset("speech", n = 32000, sr = 16000)
#'
#' # Same, but one feature vector for the whole trial.
#' scat_preset("speech", n = 32000, sr = 16000, T = "global")
#' @export
scat_preset <- function(kind = c("kinematic", "speech", "audio"), n, sr, ...) {
  kind <- match.arg(kind)
  p <- PRESETS[[kind]]
  n <- check_count(n, "n")
  if (!is.numeric(sr) || length(sr) != 1L || sr <= 0) {
    stop("sr must be a single positive number.", call. = FALSE)
  }

  # Round rather than truncate: scales are dyadic, so the nearest power of two
  # to the target is what "about 4 seconds" should mean.
  J_want <- round(log2(p$max_scale_sec * sr))
  J_max <- floor(log2(n))
  if (J_want < 1L) {
    stop("At ", sr, " Hz a ", p$max_scale_sec, " second scale is less than one ",
         "sample. This preset does not suit that sampling rate.", call. = FALSE)
  }
  J <- min(J_want, J_max)
  if (J < J_want) {
    warning(
      "The \"", kind, "\" preset wants a maximum scale of ", p$max_scale_sec,
      " s (J = ", J_want, "), but a signal of ", n, " samples only supports ",
      "J = ", J, " (", signif(2^J / sr, 3), " s). Low-frequency structure ",
      "beyond that will not be represented.",
      call. = FALSE
    )
  }

  defaults <- list(n = n, J = J, Q = p$Q, T_sec = p$T_sec, sr = sr)
  user <- list(...)
  # T and T_sec set the same thing, so let either one from the caller win.
  if ("T" %in% names(user)) {
    defaults$T_sec <- NULL
  }
  do.call(scattering_1d, utils::modifyList(defaults, user))
}


#' Make a signal the length an operator expects
#'
#' Operators are built for one signal length, but recordings rarely all have it.
#' This trims or extends a signal to `n`, which is usually the right thing to do
#' once and explicitly rather than to leave implicit in an analysis script.
#'
#' Note that the three methods are not interchangeable. `"trim"` and `"pad"`
#' preserve the time scale, so frequencies mean the same thing across trials;
#' `"resample"` does not, and will make a fast tapper and a slow tapper look
#' alike. Prefer `"trim"` or `"pad"` unless you specifically want to normalise
#' out duration.
#'
#' @param x A numeric vector, or a matrix whose columns are channels.
#' @param n Target length in samples.
#' @param how `"trim"` cuts or, if the signal is short, reflection-pads it.
#'   `"pad"` extends with `value` and never cuts, erroring if `x` is too long.
#'   `"resample"` interpolates onto `n` points, changing the time scale.
#' @param where For `"trim"`, whether to keep the `"centre"`, `"start"` or
#'   `"end"` of the signal.
#' @param value Constant used by `"pad"`.
#'
#' @return A matrix with `n` rows, or a vector if `x` was a vector.
#'
#' @examples
#' scat_fit_length(1:100, n = 64) |> length()
#' scat_fit_length(1:50, n = 64, how = "pad", value = 0) |> length()
#' @export
scat_fit_length <- function(x, n, how = c("trim", "pad", "resample"),
                            where = c("centre", "start", "end"), value = 0) {
  how <- match.arg(how)
  where <- match.arg(where)
  n <- check_count(n, "n")
  was_vector <- is.null(dim(x))
  m <- as_signal_matrix_lenient(x)
  len <- nrow(m)

  out <- if (how == "resample") {
    apply(m, 2L, function(col) {
      stats::approx(seq_along(col), col, n = n)$y
    })
  } else if (len == n) {
    m
  } else if (len > n) {
    if (how == "pad") {
      stop("x has ", len, " samples, which is more than n = ", n,
           ". Use how = \"trim\".", call. = FALSE)
    }
    from <- switch(where,
      centre = (len - n) %/% 2L + 1L,
      start = 1L,
      end = len - n + 1L
    )
    m[seq.int(from, from + n - 1L), , drop = FALSE]
  } else {
    if (how == "pad") {
      rbind(m, matrix(value, nrow = n - len, ncol = ncol(m)))
    } else {
      apply(m, 2L, function(col) {
        # Reflection keeps the local spectrum plausible where a constant would
        # inject a step.
        pad_reflect(col, 0L, n - len)
      })
    }
  }

  out <- matrix(out, nrow = n, dimnames = list(NULL, colnames(m)))
  if (was_vector) as.numeric(out) else out
}


#' Fill short gaps in a tracked signal
#'
#' Landmark tracking drops frames. Scattering has no notion of a missing sample,
#' so gaps must be dealt with before transforming. This interpolates runs of
#' missing values up to `max_gap` long and leaves anything longer alone, so that
#' a genuine tracking failure stays visible as `NA` rather than being silently
#' invented.
#'
#' @param x A numeric vector, or a matrix whose columns are channels.
#' @param max_gap Longest run of consecutive missing values to interpolate, in
#'   samples. Runs longer than this are left as `NA`.
#' @param method `"linear"` or `"spline"`.
#'
#' @return An object shaped like `x`, with short gaps filled.
#'
#' @examples
#' x <- sin(seq(0, 10, length.out = 100))
#' x[c(20, 21, 60:75)] <- NA
#' filled <- scat_fill_gaps(x, max_gap = 5)
#' sum(is.na(x))
#' sum(is.na(filled))   # the 16-sample gap is left alone
#' @export
scat_fill_gaps <- function(x, max_gap = 5, method = c("linear", "spline")) {
  method <- match.arg(method)
  max_gap <- check_count(max_gap, "max_gap")
  was_vector <- is.null(dim(x))
  m <- as_signal_matrix_lenient(x)

  filled <- apply(m, 2L, function(col) fill_gaps_one(col, max_gap, method))
  filled <- matrix(filled, nrow = nrow(m), dimnames = list(NULL, colnames(m)))
  if (was_vector) as.numeric(filled) else filled
}

#' @noRd
fill_gaps_one <- function(col, max_gap, method) {
  na <- is.na(col)
  if (!any(na) || all(na)) {
    return(col)
  }
  known <- which(!na)
  interp <- if (method == "linear") {
    stats::approx(known, col[known], xout = seq_along(col), rule = 1)$y
  } else {
    stats::spline(known, col[known], xout = seq_along(col))$y
  }
  # Only accept the interpolation inside runs that are short enough, and never
  # outside the observed range, where there is nothing to interpolate between.
  runs <- rle(na)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  ok <- rep(FALSE, length(col))
  for (i in which(runs$values & runs$lengths <= max_gap)) {
    if (starts[i] > 1L && ends[i] < length(col)) {
      ok[starts[i]:ends[i]] <- TRUE
    }
  }
  col[ok] <- interp[ok]
  col
}

#' Like as_signal_matrix() but tolerant of NA and of any length
#' @noRd
as_signal_matrix_lenient <- function(x) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = 1L)
  }
  if (!is.numeric(x)) {
    stop("x must be numeric.", call. = FALSE)
  }
  x
}
