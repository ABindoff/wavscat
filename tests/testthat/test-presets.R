test_that("presets pick scales appropriate to the sampling rate", {
  tap <- scat_preset("kinematic", n = 1800, sr = 30)
  expect_s3_class(tap, "wavscat_op")
  # A maximum scale of about 4 s at 30 Hz.
  expect_equal(tap$J, 7L)
  expect_equal(2^tap$J / 30, 4.267, tolerance = 1e-3)
  expect_equal(tap$T, 30)
  expect_equal(tap$Q, c(8L, 1L))
  # The band reaches below a plausible tapping rate and up towards Nyquist.
  xi1 <- tap$paths$xi1[tap$paths$order == 1L] * 30
  expect_lt(min(xi1), 2)
  expect_gt(max(xi1), 7)

  ddk <- scat_preset("speech", n = 32000, sr = 16000)
  expect_equal(ddk$J, 12L)
  expect_equal(ddk$T, 400)
  expect_equal(ddk$Q, c(8L, 2L))
})


test_that("presets accept overrides", {
  op <- scat_preset("kinematic", n = 1800, sr = 30, Q = 4, max_order = 1)
  expect_equal(op$Q, c(4L, 1L))
  expect_equal(op$max_order, 1L)

  # T from the caller wins over the preset's T_sec.
  g <- scat_preset("speech", n = 32000, sr = 16000, T = "global")
  expect_equal(g$average, "global")
})


test_that("a preset warns when the signal is too short for its scale", {
  # 100 samples at 30 Hz supports J = 6, but the preset wants J = 7.
  expect_warning(scat_preset("kinematic", n = 100, sr = 30, T = 8), "only supports")
})


test_that("an averaging scale wider than the signal is flagged", {
  expect_warning(scattering_1d(n = 64, J = 5, T = 32), "wider than the signal")
})


test_that("scat_fit_length trims, pads and resamples to the target", {
  x <- as.numeric(1:100)
  expect_length(scat_fit_length(x, 64), 64)
  expect_length(scat_fit_length(x, 64, how = "resample"), 64)
  expect_length(scat_fit_length(1:50, 64, how = "pad"), 64)
  expect_length(scat_fit_length(1:50, 64, how = "trim"), 64)

  # Trimming takes the requested part of the signal.
  expect_equal(scat_fit_length(x, 10, where = "start"), as.numeric(1:10))
  expect_equal(scat_fit_length(x, 10, where = "end"), as.numeric(91:100))
  expect_equal(scat_fit_length(x, 10, where = "centre"), as.numeric(46:55))

  # Padding keeps the original at the front.
  padded <- scat_fit_length(as.numeric(1:50), 64, how = "pad", value = -1)
  expect_equal(padded[1:50], as.numeric(1:50))
  expect_true(all(padded[51:64] == -1))

  expect_error(scat_fit_length(1:100, 64, how = "pad"), "more than n")
})


test_that("scat_fit_length keeps matrices and their column names", {
  m <- cbind(a = 1:100, b = 101:200)
  out <- scat_fit_length(m, 64)
  expect_equal(dim(out), c(64L, 2L))
  expect_equal(colnames(out), c("a", "b"))
})


test_that("scat_fill_gaps fills short runs and leaves long ones", {
  x <- sin(seq(0, 10, length.out = 100))
  truth <- x
  x[c(20, 21)] <- NA        # short gap, should be filled
  x[60:75] <- NA            # 16 samples, longer than max_gap
  filled <- scat_fill_gaps(x, max_gap = 5)

  expect_false(anyNA(filled[20:21]))
  expect_true(all(is.na(filled[60:75])))
  expect_lt(max(abs(filled[20:21] - truth[20:21])), 0.05)
  # Untouched samples stay exactly as they were.
  expect_equal(filled[1:19], x[1:19])
})


test_that("scat_fill_gaps leaves gaps at the ends alone", {
  x <- c(NA, NA, as.numeric(3:100), NA)
  filled <- scat_fill_gaps(x, max_gap = 5)
  expect_true(all(is.na(filled[1:2])))
  expect_true(is.na(filled[101]))
})


test_that("scat_fill_gaps handles matrices and the all-or-nothing cases", {
  m <- cbind(a = c(1, NA, 3), b = c(NA_real_, NA_real_, NA_real_))
  out <- scat_fill_gaps(m, max_gap = 3)
  expect_equal(out[, "a"], c(1, 2, 3))
  expect_true(all(is.na(out[, "b"])))
  expect_equal(colnames(out), c("a", "b"))
})


test_that("gap filling makes a tracked signal transformable", {
  set.seed(3)
  x <- sin(2 * pi * 0.05 * seq_len(1800))
  x[sample(2:1799, 40)] <- NA
  op <- scat_preset("kinematic", n = 1800, sr = 30)
  expect_error(scat_transform(op, x), "missing values")
  expect_s3_class(scat_transform(op, scat_fill_gaps(x)), "wavscat_coefs")
})
