# Plots need a graphics device, so each test opens a temporary png one.

test_that("the filter bank plot draws without error", {
  op <- scattering_1d(n = 8192, J = 8, Q = 8, sr = 16000)
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit({
    grDevices::dev.off()
    unlink(f)
  })

  expect_invisible(plot(op))
  expect_silent(plot(op, lp = TRUE))
  expect_silent(plot(op, order = 2))
  expect_error(plot(op, order = 3), "must be 1 or 2")
})


test_that("the filter bank plot works without a sampling rate", {
  op <- scattering_1d(n = 2048, J = 6, Q = 4)
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit({
    grDevices::dev.off()
    unlink(f)
  })
  expect_silent(plot(op))
})


test_that("the scalogram draws for both orders and validates its arguments", {
  t <- seq_len(8192)
  x <- sin(2 * pi * (0.002 + 0.2 * t / 8192) * t)
  S <- scattering(cbind(one = x, two = rev(x)), J = 8, Q = 8, sr = 8000)

  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit({
    grDevices::dev.off()
    unlink(f)
  })

  expect_invisible(plot(S))
  expect_silent(plot(S, order = 2))
  expect_silent(plot(S, channel = "two"))
  expect_silent(plot(scat_log(S, eps = scat_eps(S))))

  expect_error(plot(S, channel = "nope"), "channel must name or index")
  expect_error(plot(S, order = 7), "no order-7 coefficients")
})


test_that("the scalogram refuses input it cannot show", {
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit({
    grDevices::dev.off()
    unlink(f)
  })

  glob <- scattering(stats::rnorm(2048), J = 6, T = "global")
  expect_error(plot(glob), "only one time point")

  op <- scattering_1d(n = 2048, J = 6, T = 0, out_type = "list")
  lst <- scat_transform(op, stats::rnorm(2048))
  expect_error(plot(lst), "array storage")
})
