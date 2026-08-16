# Behaviour of joint time-frequency scattering, independent of Kymatio.

test_that("the frequential bank is balanced between spins", {
  op <- scattering_jtfs(n = 4096, J = 6, J_fr = 3, Q = 8)
  meta <- scat_meta(op)
  # Order two carries both spins in equal number. Order one is real-valued and
  # so only carries the non-negative ones, which is checked separately below.
  spins <- table(meta$spin[meta$order == 2L])
  expect_equal(as.integer(spins[["-1"]]), as.integer(spins[["1"]]))
  expect_gt(as.integer(spins[["0"]]), 0L)

  # Spin is the sign of the frequential centre frequency, by definition.
  expect_equal(meta$spin[meta$order > 0L], sign(meta$xi_fr[meta$order > 0L]))
})


test_that("only the non-negative spins appear at first order", {
  # S1 is real, so its negative-spin frequential responses are conjugates of
  # the positive ones and carry nothing new.
  op <- scattering_jtfs(n = 4096, J = 6, J_fr = 3, Q = 8)
  meta <- scat_meta(op)
  expect_true(all(meta$spin[meta$order == 1L] >= 0))
  expect_true(any(meta$spin[meta$order == 2L] < 0))
})


test_that("spin separates rising from falling chirps", {
  # The property that motivates joint scattering. Time scattering gives these
  # two signals nearly identical coefficients; here they land on opposite spins.
  n <- 8192
  t <- seq_len(n)
  up <- sin(2 * pi * (0.005 + 0.35 * t / n) * t)
  down <- rev(up)

  op <- scattering_jtfs(n = n, J = 8, J_fr = 3, Q = 8)
  meta <- scat_meta(op)
  o2 <- meta$order == 2L
  energy <- function(x, sp) sum(scat_transform(op, x)$coef[o2 & meta$spin == sp, , 1L]^2)

  ratio_up <- energy(up, 1) / energy(up, -1)
  ratio_down <- energy(down, 1) / energy(down, -1)

  expect_gt(ratio_up, 3)      # rising sweep favours one spin
  expect_lt(ratio_down, 1 / 3) # falling sweep favours the other
  expect_gt(ratio_up / ratio_down, 10)
})


test_that("global averaging destroys the spin asymmetry", {
  # Worth pinning down, because it is a trap. Global averaging sums over the
  # reflection-padded signal, and the mirror of a rising sweep is a falling one,
  # so the two spins end up with equal energy. This matches Kymatio, and is why
  # the documentation steers towards a local T.
  n <- 8192
  t <- seq_len(n)
  up <- sin(2 * pi * (0.005 + 0.35 * t / n) * t)

  op <- scattering_jtfs(n = n, J = 8, J_fr = 3, Q = 8, T = "global")
  meta <- scat_meta(op)
  o2 <- meta$order == 2L
  S <- scat_transform(op, up)$coef[, 1L, 1L]
  ratio <- sum(S[o2 & meta$spin == 1]^2) / sum(S[o2 & meta$spin == -1]^2)
  expect_lt(abs(ratio - 1), 0.05)
})


test_that("joint scattering is homogeneous and non-negative", {
  set.seed(21)
  op <- scattering_jtfs(n = 2048, J = 5, J_fr = 2, Q = 4)
  x <- stats::rnorm(2048)
  Sx <- scat_transform(op, x)$coef
  Sax <- scat_transform(op, 4.5 * x)$coef
  expect_lt(max(abs(Sax - 4.5 * Sx)) / max(abs(4.5 * Sx)), 1e-12)

  meta <- scat_meta(op)
  expect_gte(min(Sx[meta$order >= 1L, , 1L]), 0)
})


test_that("the two output formats agree where they overlap", {
  set.seed(22)
  x <- stats::rnorm(4096)
  args <- list(n = 4096, J = 6, J_fr = 3, Q = 8)

  time_fmt <- scat_transform(do.call(scattering_jtfs, c(args, format = "time")), x)
  joint <- scat_transform(
    do.call(scattering_jtfs, c(args, list(format = "joint", out_type = "list"))), x
  )
  mj <- joint$meta
  mt <- time_fmt$meta

  # Each time-format path is one band of the corresponding joint path.
  row <- which(mt$order == 1L & mt$n_fr == 2L)
  jrow <- which(mj$order == 1L & mj$n_fr == 2L)
  band_of <- joint$coef[[jrow]][, , 1L]
  for (k in seq_along(row)) {
    expect_equal(time_fmt$coef[row[k], , 1L], band_of[k, ], tolerance = 1e-12)
  }
})


test_that("joint output keeps a frequency axis and drops order zero", {
  op <- scattering_jtfs(n = 4096, J = 6, J_fr = 3, Q = 8, format = "joint")
  S <- scat_transform(op, stats::rnorm(4096))
  expect_length(dim(S$coef), 4L)
  expect_equal(nrow(S$meta), dim(S$coef)[1L])
  # Order zero has no frequency axis, so an array of joint paths cannot hold it.
  expect_false(any(S$meta$order == 0L))

  time_fmt <- scattering_jtfs(n = 4096, J = 6, J_fr = 3, Q = 8, format = "time")
  St <- scat_transform(time_fmt, stats::rnorm(4096))
  expect_length(dim(St$coef), 3L)
  expect_true(any(St$meta$order == 0L))
})


test_that("frequential averaging controls the number of paths", {
  count <- function(...) {
    nrow(scat_meta(scattering_jtfs(n = 4096, J = 6, J_fr = 3, Q = 8, ...)))
  }
  # Coarser frequential averaging pools more bands, so fewer paths survive.
  expect_gt(count(F = 0), count(F = 4))
  expect_gt(count(F = 4), count(F = "global"))
})


test_that("joint operators validate their arguments", {
  expect_error(
    scattering_jtfs(n = 4096, J = 6, J_fr = 3, F = 0, format = "joint"),
    "needs frequential averaging"
  )
  expect_error(scattering_jtfs(n = 4096, J = 6, J_fr = 0), "positive whole number")
  expect_error(scattering_jtfs(n = 4096, J = 6, J_fr = 9), "only 56 first-order bands")
  expect_error(
    scattering_jtfs(n = 4096, J = 6, J_fr = 3, stride_fr = 3),
    "power of two"
  )
})


test_that("joint metadata reports frequential rate per octave", {
  op <- scattering_jtfs(n = 4096, J = 6, J_fr = 3, Q = 8, sr = 16000)
  meta <- scat_meta(op)
  expect_equal(meta$xi_fr_octave, meta$xi_fr * 8)
  expect_true(all(is.na(meta$xi_fr) | abs(meta$xi_fr) <= 0.5))
})


test_that("joint operators print and size themselves", {
  op <- scattering_jtfs(n = 4096, J = 6, J_fr = 3, Q = 8, sr = 16000)
  expect_output(print(op), "joint time-frequency")
  expect_output(print(op), "J_fr")
  expect_equal(scat_output_size(op), nrow(scat_meta(op)))
})
