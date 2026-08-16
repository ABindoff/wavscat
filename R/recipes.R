# A recipes step, so that scattering can live inside a resampled workflow.

#' Scattering features in a recipe
#'
#' `step_scattering()` creates a specification of a recipe step that replaces a
#' block of columns holding one signal with its scattering features.
#'
#' The data is expected wide: one row per recording, and the selected columns
#' are the time samples of that recording, in order. Select them with something
#' like `starts_with("t")` or `all_numeric_predictors()`.
#'
#' @section Why put it in a recipe:
#'
#' The compression scale `eps` has to come from somewhere, and choosing it from
#' the whole dataset leaks test information into training. Inside a recipe,
#' `prep()` sees only the analysis half of each resample, so `eps` is learned
#' there and applied unchanged when baking the assessment half. Doing the
#' transform by hand before splitting quietly gets this wrong.
#'
#' @param recipe A recipe object.
#' @param ... One or more selector functions choosing the columns that hold the
#'   signal. See [recipes::selections()].
#' @param role Role assigned to the new columns.
#' @param trained Whether the step has been prepared.
#' @param J,Q,T,T_sec,max_order,sr Passed to [scattering_1d()]. `J` defaults to
#'   the largest scale the signal length supports, minus one.
#' @param renorm Whether to divide order-2 coefficients by their order-1 parent,
#'   as in [scat_renorm()].
#' @param log Whether to log-compress, as in [scat_log()].
#' @param eps Compression scale. Leave `NULL` to learn it from the training data
#'   during `prep()`, which is the point of using a recipe.
#' @param summary How to reduce the time axis, passed to [scat_features()].
#' @param prefix Prepended to the new column names.
#' @param n_eps Number of training rows sampled to estimate `eps`. Estimating it
#'   from every row is usually a waste; a few dozen fix the order of magnitude,
#'   which is all `eps` needs.
#' @param operator,columns,eps_value Populated during `prep()`; not user-facing.
#' @param keep_original_cols Whether to keep the input columns as well.
#' @param skip Whether to skip this step when baking new data. Leave `FALSE`.
#' @param id A unique identifier for the step.
#'
#' @return An updated recipe with the new step added.
#'
#' @examplesIf requireNamespace("recipes", quietly = TRUE)
#' set.seed(1)
#' # Twenty trials of 512 samples each, wide: one row per trial.
#' trials <- as.data.frame(matrix(rnorm(20 * 512), nrow = 20))
#' trials$group <- rep(c("a", "b"), each = 10)
#'
#' rec <- recipes::recipe(group ~ ., data = trials) |>
#'   step_scattering(recipes::all_numeric_predictors(), J = 5, Q = 4)
#'
#' prepped <- recipes::prep(rec, training = trials)
#' out <- recipes::bake(prepped, new_data = NULL)
#' dim(out)
#' names(out)[1:4]
#' @export
step_scattering <- function(recipe, ..., role = "predictor", trained = FALSE,
                            J = NULL, Q = 8, T = NULL, T_sec = NULL,
                            max_order = 2, sr = NULL,
                            renorm = TRUE, log = TRUE, eps = NULL,
                            summary = "mean", prefix = "scat_", n_eps = 50,
                            operator = NULL, columns = NULL, eps_value = NULL,
                            keep_original_cols = FALSE, skip = FALSE,
                            id = recipes::rand_id("scattering")) {
  recipes::add_step(
    recipe,
    step_scattering_new(
      terms = rlang::enquos(...), role = role, trained = trained,
      J = J, Q = Q, T = T, T_sec = T_sec, max_order = max_order, sr = sr,
      renorm = renorm, log = log, eps = eps, summary = summary,
      prefix = prefix, n_eps = n_eps,
      operator = operator, columns = columns, eps_value = eps_value,
      keep_original_cols = keep_original_cols, skip = skip, id = id
    )
  )
}

#' @noRd
step_scattering_new <- function(terms, role, trained, J, Q, T, T_sec, max_order,
                                sr, renorm, log, eps, summary, prefix, n_eps,
                                operator, columns, eps_value,
                                keep_original_cols, skip, id) {
  recipes::step(
    subclass = "scattering",
    terms = terms, role = role, trained = trained,
    J = J, Q = Q, T = T, T_sec = T_sec, max_order = max_order, sr = sr,
    renorm = renorm, log = log, eps = eps, summary = summary,
    prefix = prefix, n_eps = n_eps,
    operator = operator, columns = columns, eps_value = eps_value,
    keep_original_cols = keep_original_cols, skip = skip, id = id
  )
}

#' @noRd
prep.step_scattering <- function(x, training, info = NULL, ...) {
  col_names <- recipes::recipes_eval_select(x$terms, training, info)
  recipes::check_type(training[, col_names], types = c("double", "integer"))

  n <- length(col_names)
  if (n < 8L) {
    stop("step_scattering() selected ", n, " column", if (n != 1L) "s" else "",
         ". The selected columns should be the time samples of one signal, so ",
         "there need to be many of them; check the selector.", call. = FALSE)
  }

  J <- x$J %||% max(floor(log2(n)) - 1L, 1L)
  op <- scattering_1d(
    n = n, J = J, Q = x$Q, T = x$T, T_sec = x$T_sec,
    max_order = x$max_order, sr = x$sr
  )

  eps_value <- x$eps
  if (isTRUE(x$log) && is.null(eps_value)) {
    rows <- unique(round(seq(1, nrow(training), length.out = min(x$n_eps, nrow(training)))))
    sample_coefs <- scat_rows(op, training[rows, col_names, drop = FALSE],
                              renorm = x$renorm)
    eps_value <- scat_eps(sample_coefs)
  }

  step_scattering_new(
    terms = x$terms, role = x$role, trained = TRUE,
    J = J, Q = x$Q, T = x$T, T_sec = x$T_sec, max_order = x$max_order,
    sr = x$sr, renorm = x$renorm, log = x$log, eps = x$eps,
    summary = x$summary, prefix = x$prefix, n_eps = x$n_eps,
    operator = op, columns = col_names, eps_value = eps_value,
    keep_original_cols = x$keep_original_cols, skip = x$skip, id = x$id
  )
}

#' @noRd
bake.step_scattering <- function(object, new_data, ...) {
  coefs <- scat_rows(object$operator, new_data[, object$columns, drop = FALSE],
                     renorm = object$renorm)
  if (isTRUE(object$log)) {
    coefs <- scat_log(coefs, eps = object$eps_value)
  }
  feats <- scat_features(coefs, summary = object$summary, prefix = object$prefix)
  feats$channel <- NULL

  if (!object$keep_original_cols) {
    new_data <- new_data[, !(names(new_data) %in% object$columns), drop = FALSE]
  }
  tibble::as_tibble(cbind(feats, new_data))
}

#' @export
print.step_scattering <- function(x, width = max(20, options()$width - 35), ...) {
  title <- if (x$trained) {
    sprintf("Wavelet scattering (J = %d, Q = %s) on ", x$J,
            paste(x$operator$Q, collapse = "/"))
  } else {
    "Wavelet scattering on "
  }
  recipes::print_step(x$columns, x$terms, x$trained, title, width)
  invisible(x)
}

#' @noRd
tidy.step_scattering <- function(x, ...) {
  if (x$trained) {
    res <- tibble::tibble(
      terms = x$columns,
      J = x$J,
      eps = x$eps_value %||% NA_real_
    )
  } else {
    res <- tibble::tibble(
      terms = recipes::sel2char(x$terms),
      J = x$J %||% NA_integer_,
      eps = NA_real_
    )
  }
  res$id <- x$id
  res
}

#' @noRd
required_pkgs.step_scattering <- function(x, ...) {
  "wavscat"
}

#' Transform a wide block of signals, one row per recording
#'
#' Rows become channels of a single call, which is both simpler and faster than
#' looping: the padded spectrum machinery is set up once per path rather than
#' once per row.
#'
#' @noRd
scat_rows <- function(op, data, renorm = FALSE) {
  m <- t(as.matrix(data))
  colnames(m) <- rownames(data) %||% as.character(seq_len(ncol(m)))
  coefs <- scat_transform(op, m)
  if (isTRUE(renorm)) {
    coefs <- scat_renorm(coefs)
  }
  coefs
}
