skip_if_not_installed("recipes")

make_trials <- function(n_trial = 20, n_time = 512, seed = 1) {
  set.seed(seed)
  d <- as.data.frame(matrix(stats::rnorm(n_trial * n_time), nrow = n_trial))
  names(d) <- paste0("t", seq_len(n_time))
  d$group <- factor(rep(c("a", "b"), length.out = n_trial))
  d
}


test_that("the step prepares and bakes to scattering features", {
  trials <- make_trials()
  rec <- recipes::recipe(group ~ ., data = trials) |>
    step_scattering(recipes::all_numeric_predictors(), J = 5, Q = 4)

  prepped <- recipes::prep(rec, training = trials)
  out <- recipes::bake(prepped, new_data = NULL)

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), nrow(trials))
  # The 512 input columns are gone, replaced by features plus the outcome.
  expect_false(any(startsWith(names(out), "t1")))
  expect_true("group" %in% names(out))
  expect_true(all(startsWith(setdiff(names(out), "group"), "scat_")))
  expect_false(anyNA(out))
})


test_that("features match calling the functions directly", {
  trials <- make_trials()
  cols <- paste0("t", seq_len(512))
  rec <- recipes::recipe(group ~ ., data = trials) |>
    step_scattering(recipes::all_numeric_predictors(), J = 5, Q = 4, eps = 1e-4)
  out <- recipes::bake(recipes::prep(rec, training = trials), new_data = NULL)

  op <- scattering_1d(n = 512, J = 5, Q = 4)
  S <- scat_log(scat_renorm(scat_transform(op, t(as.matrix(trials[, cols])))),
                eps = 1e-4)
  direct <- scat_features(S, prefix = "scat_")
  direct$channel <- NULL

  expect_equal(as.data.frame(out[, names(direct)]), as.data.frame(direct))
})


test_that("eps is learned from training data and reused when baking", {
  trials <- make_trials(n_trial = 30)
  rec <- recipes::recipe(group ~ ., data = trials) |>
    step_scattering(recipes::all_numeric_predictors(), J = 5, Q = 4)

  # Prepare on one half; the eps must not change when the other half is baked.
  prepped <- recipes::prep(rec, training = trials[1:15, ])
  learned <- recipes::tidy(prepped, number = 1)$eps[1]
  expect_true(is.finite(learned) && learned > 0)

  baked_a <- recipes::bake(prepped, new_data = trials[1:15, ])
  baked_b <- recipes::bake(prepped, new_data = trials[16:30, ])
  expect_equal(ncol(baked_a), ncol(baked_b))
  expect_equal(names(baked_a), names(baked_b))

  # Baking a single row gives the same values as baking it inside a block,
  # which is what guarantees no information crosses between rows.
  one <- recipes::bake(prepped, new_data = trials[16, ])
  expect_equal(as.numeric(one[1, ]), as.numeric(baked_b[1, ]))
})


test_that("an explicit eps overrides the learned one", {
  trials <- make_trials()
  rec <- recipes::recipe(group ~ ., data = trials) |>
    step_scattering(recipes::all_numeric_predictors(), J = 5, Q = 4, eps = 0.25)
  prepped <- recipes::prep(rec, training = trials)
  expect_equal(recipes::tidy(prepped, number = 1)$eps[1], 0.25)
})


test_that("log and renorm can be switched off", {
  trials <- make_trials()
  rec <- recipes::recipe(group ~ ., data = trials) |>
    step_scattering(recipes::all_numeric_predictors(), J = 5, Q = 4,
                    log = FALSE, renorm = FALSE)
  prepped <- recipes::prep(rec, training = trials)
  out <- recipes::bake(prepped, new_data = NULL)

  op <- scattering_1d(n = 512, J = 5, Q = 4)
  S <- scat_transform(op, t(as.matrix(trials[, paste0("t", 1:512)])))
  direct <- scat_features(S, prefix = "scat_")
  direct$channel <- NULL
  expect_equal(as.data.frame(out[, names(direct)]), as.data.frame(direct))
  expect_true(is.na(recipes::tidy(prepped, number = 1)$eps[1]))
})


test_that("the original columns can be kept", {
  trials <- make_trials(n_time = 256)
  rec <- recipes::recipe(group ~ ., data = trials) |>
    step_scattering(recipes::all_numeric_predictors(), J = 4, Q = 4,
                    keep_original_cols = TRUE)
  out <- recipes::bake(recipes::prep(rec, training = trials), new_data = NULL)
  expect_true("t1" %in% names(out))
})


test_that("the step prints and tidies before and after prep", {
  trials <- make_trials(n_time = 256)
  rec <- recipes::recipe(group ~ ., data = trials) |>
    step_scattering(recipes::all_numeric_predictors(), J = 4, Q = 4)
  # recipes prints through cli, so capture both streams.
  shown <- function(x) {
    paste(c(utils::capture.output(print(x)),
            utils::capture.output(print(x), type = "message")), collapse = "\n")
  }
  expect_match(shown(rec), "Wavelet scattering")

  untrained <- recipes::tidy(rec, number = 1)
  expect_s3_class(untrained, "tbl_df")
  expect_true(all(is.na(untrained$eps)))

  prepped <- recipes::prep(rec, training = trials)
  expect_match(shown(prepped), "Wavelet scattering")
  expect_true("wavscat" %in% recipes::required_pkgs(prepped))
})


test_that("selecting too few columns is caught", {
  d <- data.frame(a = 1:10, b = 11:20, group = factor(rep(1:2, 5)))
  rec <- recipes::recipe(group ~ ., data = d) |>
    step_scattering(recipes::all_numeric_predictors())
  expect_error(recipes::prep(rec, training = d), "check the selector")
})
