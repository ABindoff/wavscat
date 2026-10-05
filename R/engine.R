# The optional Rust engine.
#
# wavscat computes everything in base R. If the companion package
# 'wavscatengine' is installed, it computes the coefficients instead, with
# the same Rust core used from Python and in the browser, so that results are
# bit-for-bit identical across languages and platforms. The two engines agree
# to about 1e-14 either way; only the Rust one promises identity.
#
# Choose with options(wavscat.engine = ...): "auto" (the default) uses Rust
# when it is installed, "rust" requires it, "r" never uses it.

#' Which engine computes the coefficients
#'
#' wavscat computes scattering coefficients in base R. When the companion
#' package 'wavscatengine' is installed, it uses that instead: the same Rust
#' code that runs in Python and in the browser, so features are bit-for-bit
#' identical across languages and platforms. The two engines agree to about
#' 1e-14 in any case; only the Rust one guarantees identity.
#'
#' Choose with `options(wavscat.engine = "auto")` (the default: Rust when it
#' is installed), `"rust"` (require it) or `"r"` (never use it). Each set of
#' coefficients records the engine that produced it in `x$spec$engine`.
#'
#' @return `"rust"` or `"r"`, the engine the next transform will use.
#' @examples
#' scat_engine()
#' @export
scat_engine <- function() {
  choice <- getOption("wavscat.engine", "auto")
  if (!choice %in% c("auto", "rust", "r")) {
    stop("options(wavscat.engine) must be \"auto\", \"rust\" or \"r\".", call. = FALSE)
  }
  available <- requireNamespace("wavscatengine", quietly = TRUE)
  if (choice == "rust" && !available) {
    stop("options(wavscat.engine = \"rust\") needs the 'wavscatengine' package, ",
         "which is not installed.", call. = FALSE)
  }
  if (choice == "r" || !available) "r" else "rust"
}

#' Resolved parameters of an operator, as the engine expects them
#'
#' Positions are fixed by the engine's C interface; every value is the one
#' this operator resolved, so both engines build the same operator.
#' @noRd
engine_params <- function(op) {
  t_kind <- if (identical(op$average, "global")) 2 else 1
  t_value <- if (isFALSE(op$average)) 0 else op$T
  stride <- if (identical(op$average, "local")) 2^op$log2_stride else NaN
  if (!identical(op$type, "jtfs")) {
    return(c(op$n, op$J, op$Q[1L], op$Q[2L], t_kind, t_value, op$max_order, stride))
  }
  f_kind <- if (identical(op$average_fr, "global")) 2 else 1
  f_value <- if (isFALSE(op$average_fr)) 0 else op$F
  stride_fr <- if (identical(op$average_fr, "local")) 2^op$log2_stride_fr else NaN
  c(op$n, op$J, op$Q[1L], op$Q[2L], t_kind, t_value, stride, op$J_fr, op$Q_fr,
    f_kind, f_value, stride_fr, if (identical(op$format, "joint")) 1 else 0,
    if (identical(op$out_type, "list")) 1 else 0)
}

#' Run the Rust engine over every channel
#'
#' @return A list with `per_channel`, in the shape `scat_run()` gives, and for
#'   joint scattering `s1`, one first-order matrix per channel.
#' @noRd
engine_run <- function(op, x) {
  params <- engine_params(op)
  jtfs <- identical(op$type, "jtfs")
  r <- if (jtfs) {
    wavscatengine::engine_jtfs(params, x)
  } else {
    wavscatengine::engine_scattering_1d(params, x)
  }
  if (!identical(r$labels, scat_meta(op)$path)) {
    stop("The Rust engine produced paths in a different order from this ",
         "operator. Please report this as a bug.", call. = FALSE)
  }
  orders <- op$paths$order
  keep_matrix <- jtfs && identical(op$format, "joint")
  per_channel <- lapply(r$coef, function(paths) {
    Map(function(coef, order) {
      # Joint paths keep their [band, time] shape, except order zero, which
      # has no frequency axis; everything else is a plain series.
      if (!keep_matrix || order == 0L) coef <- as.vector(coef)
      list(coef = coef)
    }, paths, orders)
  })
  list(per_channel = per_channel, s1 = r$s1)
}

#' Version of the engine's numerics, or NA for the R engine
#' @noRd
engine_numerics <- function(engine) {
  if (identical(engine, "rust")) wavscatengine::engine_version()$numerics else NA_character_
}
