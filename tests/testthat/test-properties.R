# Mathematical properties of the transform.
#
# These are independent of Kymatio: a faithful transcription of a misread
# formula would still match the reference. They check that the object really is
# a scattering transform.

CONFIGS <- list(
  list(N = 2048, J = 6, Q = 8, T = 64),
  list(N = 4096, J = 8, Q = 8, T = 256),
  list(N = 2048, J = 6, Q = 1, T = 64),
  list(N = 4096, J = 8, Q = 16, T = 256),
  list(N = 1024, J = 5, Q = 4, T = 32),
  list(N = 8192, J = 10, Q = 12, T = 1024)
)


test_that("wavelets have exactly zero mean and the low-pass has unit gain", {
  for (cfg in CONFIGS) {
    fb <- scattering_filter_factory(cfg$N, cfg$J, c(cfg$Q, 1), cfg$T,
                                    5, sqrt(0.5), 0.1)
    tag <- sprintf("N=%d J=%d Q=%d", cfg$N, cfg$J, cfg$Q)

    # psi_hat(0) == 0 is what makes these wavelets rather than mere band-passes.
    for (bank in list(fb$psi1, fb$psi2)) {
      dc <- vapply(bank, function(p) p$levels[[1L]][1L], 0)
      expect_lt(max(abs(dc)), 1e-15, label = paste("wavelet DC gain,", tag))
    }
    # phi is non-negative in time with unit L1 norm, so phi_hat(0) == 1 exactly.
    expect_equal(fb$phi$levels[[1L]][1L], 1, tolerance = 1e-12,
                 label = paste("low-pass DC gain,", tag))
  }
})


test_that("the filter bank is a frame with upper bound one", {
  # The Littlewood-Paley sum |phi|^2 + (1/2) sum_lambda (|psi(w)|^2 + |psi(-w)|^2)
  # bounds the operator norm of one scattering layer. Staying at or below one is
  # exactly what makes the transform non-expansive; a bank with gaps would dip,
  # and one with too much overlap would exceed it.
  for (cfg in CONFIGS) {
    fb <- scattering_filter_factory(cfg$N, cfg$J, c(cfg$Q, 1), cfg$T,
                                    5, sqrt(0.5), 0.1)
    lp <- fb$phi$levels[[1L]]^2
    for (p in fb$psi1) {
      h <- p$levels[[1L]]
      lp <- lp + 0.5 * (h^2 + rev(c(h[1L], h[-1L]))^2)
    }
    tag <- sprintf("N=%d J=%d Q=%d", cfg$N, cfg$J, cfg$Q)
    expect_lt(max(lp), 1.001, label = paste("Littlewood-Paley upper bound,", tag))
    expect_gt(max(lp), 0.999, label = paste("Littlewood-Paley peak,", tag))
  }
})


test_that("scattering is homogeneous of degree one", {
  # Modulus is homogeneous and every other step is linear, so S(a x) = a S(x)
  # exactly for a > 0. Any indexing or normalisation slip breaks this.
  set.seed(1)
  op <- scattering_1d(n = 2048, J = 6, Q = 8)
  x <- stats::rnorm(2048)
  Sx <- scat_transform(op, x)$coef
  for (a in c(0.001, 3.7, 1000)) {
    Sax <- scat_transform(op, a * x)$coef
    expect_lt(max(abs(Sax - a * Sx)) / max(abs(a * Sx)), 1e-12,
              label = sprintf("homogeneity at a = %g", a))
  }
})


test_that("coefficients of order one and above are non-negative", {
  # They are moduli averaged against a positive kernel. Order zero is a local
  # mean and may legitimately be negative.
  set.seed(2)
  op <- scattering_1d(n = 2048, J = 6, Q = 8)
  S <- scat_transform(op, stats::rnorm(2048))
  meta <- scat_meta(op)
  expect_gte(min(S$coef[meta$order >= 1L, , 1L]), 0)
})


test_that("scattering is non-expansive", {
  set.seed(3)
  op <- scattering_1d(n = 2048, J = 6, Q = 8)
  for (i in 1:5) {
    x <- stats::rnorm(2048)
    y <- x + stats::rnorm(2048, sd = runif(1, 0.01, 2))
    Sx <- scat_transform(op, x)$coef
    Sy <- scat_transform(op, y)$coef
    expect_lte(sqrt(sum((Sx - Sy)^2)), sqrt(sum((x - y)^2)))
  }
})


test_that("energy decreases with scattering order", {
  set.seed(4)
  op <- scattering_1d(n = 2048, J = 6, Q = 8)
  for (x in list(stats::rnorm(2048),
                 sin(2 * pi * (0.005 + 0.4 * seq_len(2048) / 2048) * seq_len(2048)),
                 as.numeric(stats::filter(stats::rnorm(2048), 0.9, "recursive")))) {
    e <- scat_energy(scat_transform(op, x))
    expect_equal(e$order, 0:2)
    expect_gt(e$energy[e$order == 1L], e$energy[e$order == 2L])
  }
})


test_that("sensitivity to translation grows linearly in the shift over T", {
  # The defining property: ||S x - S(x shifted by d)|| is of order d / T for
  # shifts well inside the averaging window.
  n <- 2048
  T <- 64
  op <- scattering_1d(n = n, J = 6, Q = 8, T = T)
  tt <- seq_len(n)
  envelope <- function(d) {
    sin(2 * pi * 0.1 * (tt - d)) * exp(-((tt - d - n / 2) / 300)^2)
  }
  S0 <- scat_transform(op, envelope(0))$coef
  ref <- sqrt(sum(S0^2))

  shifts <- c(1, 2, 4, 8, 16, 32)
  rel <- vapply(shifts, function(d) {
    sqrt(sum((scat_transform(op, envelope(d))$coef - S0)^2)) / ref
  }, 0)

  # Small shifts barely move the representation at all.
  expect_lt(rel[1L], 0.01)
  # Monotone in the shift.
  expect_true(all(diff(rel) > 0))
  # Linear: doubling the shift doubles the change, to within 5 per cent.
  ratios <- rel[-1L] / rel[-length(rel)]
  expect_true(all(abs(ratios - 2) < 0.1))
})


test_that("a pure tone excites one first-order band and nothing at order two", {
  op <- scattering_1d(n = 4096, J = 8, Q = 8, T = "global")
  meta <- scat_meta(op)
  xi_true <- 0.08
  S <- scat_transform(op, sin(2 * pi * xi_true * seq_len(4096)))$coef[, 1L, 1L]

  o1 <- which(meta$order == 1L)
  peak <- o1[which.max(S[o1])]
  # The nearest filter is within half a step of the bank spacing, 2^(1/Q).
  expect_lt(abs(log2(meta$xi1[peak] / xi_true)), 1 / (2 * 8))

  # A tone has a constant envelope, so there is nothing for order two to see.
  o2 <- which(meta$order == 2L)
  expect_lt(sum(S[o2]^2) / sum(S[o1]^2), 1e-3)
})


test_that("amplitude modulation is what order two responds to", {
  op <- scattering_1d(n = 4096, J = 8, Q = 8, T = "global")
  meta <- scat_meta(op)
  o1 <- which(meta$order == 1L)
  o2 <- which(meta$order == 2L)
  carrier <- sin(2 * pi * 0.15 * seq_len(4096))

  plain <- scat_transform(op, carrier)$coef[, 1L, 1L]
  modulated <- scat_transform(
    op, carrier * (1 + 0.8 * sin(2 * pi * 0.004 * seq_len(4096)))
  )$coef[, 1L, 1L]

  ratio <- function(s) sum(s[o2]^2) / sum(s[o1]^2)
  # Modulating the carrier moves energy into order two by orders of magnitude;
  # this is the structure a first-order representation cannot see.
  expect_gt(ratio(modulated) / ratio(plain), 50)

  # The responding path should sit at the carrier frequency.
  peak <- o2[which.max(modulated[o2])]
  expect_lt(abs(log2(meta$xi1[peak] / 0.15)), 1 / (2 * 8))
})
