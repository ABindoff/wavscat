test_that("operators validate their arguments", {
  expect_error(scattering_1d(n = 100, J = 10), "exceeds the signal length")
  expect_error(scattering_1d(n = 1024, J = 5, Q = 0), "at least 1")
  expect_error(scattering_1d(n = 1024, J = 5, Q = c(1, 2, 3)), "one or two")
  expect_error(scattering_1d(n = 1024, J = 5, max_order = 3), "must be 1 or 2")
  expect_error(scattering_1d(n = 1024, J = 5, T = 2048), "exceeds the signal length")
  expect_error(scattering_1d(n = 1024, J = 5, T = 0.5), "at least 1")
  expect_error(scattering_1d(n = 1024, J = 5, T = 0), "out_type must be")
  expect_error(scattering_1d(n = 1024, J = 5, stride = 3), "power of two")
  expect_error(scattering_1d(n = 1024, J = 5, stride = 128), "coarser")
  expect_error(scattering_1d(n = 1024, J = 5, T_sec = 1), "needs sr")
  expect_error(scattering_1d(n = 1024, J = 5, sr = -1), "positive number")
})


test_that("T_sec converts seconds to samples", {
  op <- scattering_1d(n = 1800, J = 8, T_sec = 1, sr = 30)
  expect_equal(op$T, 30)
  expect_equal(op$log2_stride, 4L)
})


test_that("transform input is validated", {
  op <- scattering_1d(n = 1024, J = 5)
  expect_error(scat_transform(op, stats::rnorm(999)), "expects signals of length")
  expect_error(scat_transform(op, c(stats::rnorm(1023), NA)), "missing values")
  expect_error(scat_transform(op, as.character(seq_len(1024))), "must be numeric")
  expect_error(scat_transform("not an operator", stats::rnorm(10)), "must be a scattering operator")
})


test_that("channels are transformed independently and labelled", {
  op <- scattering_1d(n = 1024, J = 5, Q = 4)
  set.seed(7)
  x <- cbind(thumb = stats::rnorm(1024), index = stats::rnorm(1024))

  S <- scat_transform(op, x)
  expect_equal(dim(S$coef)[3L], 2L)
  expect_equal(dimnames(S$coef)$channel, c("thumb", "index"))

  for (j in 1:2) {
    solo <- scat_transform(op, x[, j])
    expect_equal(S$coef[, , j], solo$coef[, , 1L])
  }
})


test_that("a bare vector gets a default channel name", {
  S <- scattering(stats::rnorm(1024), J = 5)
  expect_equal(dimnames(S$coef)$channel, "signal")
  expect_equal(dim(S$coef)[3L], 1L)
})


test_that("data frames are accepted as multichannel input", {
  op <- scattering_1d(n = 512, J = 4, Q = 2)
  d <- data.frame(a = stats::rnorm(512), b = stats::rnorm(512))
  S <- scat_transform(op, d)
  expect_equal(dimnames(S$coef)$channel, c("a", "b"))
})


test_that("path metadata lines up with the coefficient rows", {
  op <- scattering_1d(n = 2048, J = 6, Q = 8, sr = 16000)
  meta <- scat_meta(op)
  S <- scat_transform(op, stats::rnorm(2048))

  expect_equal(nrow(meta), dim(S$coef)[1L])
  expect_equal(meta$path, dimnames(S$coef)$path)
  expect_equal(nrow(meta), scat_output_size(op))
  expect_equal(sum(scat_output_size(op, by_order = TRUE)), scat_output_size(op))

  # Order zero exactly once, orders sorted, and Hz reported when sr is known.
  expect_equal(sum(meta$order == 0L), 1L)
  expect_false(is.unsorted(meta$order))
  expect_equal(meta$xi1_hz, meta$xi1 * 16000)

  # Second-order paths always sit below their parent in frequency.
  o2 <- meta$order == 2L
  expect_true(all(meta$xi2[o2] < meta$xi1[o2]))
})


test_that("max_order = 1 drops the second order entirely", {
  op <- scattering_1d(n = 1024, J = 5, Q = 4, max_order = 1)
  meta <- scat_meta(op)
  expect_equal(max(meta$order), 1L)
  expect_true(all(is.na(meta$n2)))
})


test_that("global averaging yields one coefficient per path", {
  op <- scattering_1d(n = 1024, J = 5, Q = 4, T = "global")
  S <- scat_transform(op, stats::rnorm(1024))
  expect_equal(dim(S$coef)[2L], 1L)
  expect_true(is.na(S$spec$sr_out))
})


test_that("the output sampling rate is the input rate over the stride", {
  op <- scattering_1d(n = 4096, J = 6, Q = 8, sr = 1000)
  S <- scat_transform(op, stats::rnorm(4096))
  expect_equal(S$spec$sr_out, 1000 / 2^op$log2_stride)
})


test_that("a finer stride oversamples the output in time", {
  base <- scattering_1d(n = 2048, J = 6, Q = 8)
  fine <- scattering_1d(n = 2048, J = 6, Q = 8, stride = 2^(base$log2_stride - 2))
  Sb <- scat_transform(base, stats::rnorm(2048))
  Sf <- scat_transform(fine, stats::rnorm(2048))
  expect_equal(dim(Sf$coef)[1L], dim(Sb$coef)[1L])
  expect_gt(dim(Sf$coef)[2L], dim(Sb$coef)[2L])
})


test_that("scattering() sizes an operator to the signal", {
  S <- scattering(stats::rnorm(3000), Q = 4)
  expect_s3_class(S, "wavscat_coefs")
  expect_equal(S$spec$n, 3000L)
})


test_that("print methods run and return their input invisibly", {
  op <- scattering_1d(n = 2048, J = 6, Q = 8, sr = 16000)
  S <- scat_transform(op, stats::rnorm(2048))
  expect_output(print(op), "wavscat_op")
  expect_output(print(op), "band-pass")
  expect_output(print(S), "wavscat_coefs")
  expect_invisible(print(op))
})


test_that("reflection padding mirrors without repeating the edge", {
  expect_equal(pad_reflect(1:5, 2, 2), c(3, 2, 1, 2, 3, 4, 5, 4, 3))
  expect_equal(pad_reflect(1:5, 0, 0), 1:5)
  expect_error(pad_reflect(1:5, 5, 0), "shorter than the signal")
})


test_that("periodising a spectrum subsamples the signal", {
  set.seed(11)
  x <- stats::rnorm(64)
  for (k in c(2, 4, 8)) {
    direct <- x[seq(1, 64, by = k)]
    viafft <- Re(bk_ifft(bk_periodize(bk_fft(x), k)))
    expect_equal(viafft, direct, tolerance = 1e-12)
  }
})
