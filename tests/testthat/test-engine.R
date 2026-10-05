skip_if_not_installed("wavscatengine")

with_engine <- function(engine, code) {
  old <- options(wavscat.engine = engine)
  on.exit(options(old))
  force(code)
}

# Two channels of an irregular signal, long enough for every operator below.
signal <- function(n) {
  set.seed(11)
  cbind(a = as.numeric(arima.sim(list(ar = 0.8), n)),
        b = sin(seq_len(n) / 5) + rnorm(n, sd = 0.3))
}

operators <- list(
  "1d default" = function() scattering_1d(n = 2048, J = 6, Q = 8),
  "1d global" = function() scattering_1d(n = 2048, J = 6, Q = c(8, 2), T = "global"),
  "1d unaveraged" = function() scattering_1d(n = 1500, J = 5, Q = 4, T = 0, out_type = "list"),
  "1d strided" = function() scattering_1d(n = 2048, J = 7, Q = 8, stride = 16),
  "1d order one" = function() scattering_1d(n = 1800, J = 8, Q = 12, max_order = 1),
  "1d T_sec" = function() scattering_1d(n = 900, J = 7, Q = c(8, 1), T_sec = 4, sr = 30),
  "jtfs time" = function() scattering_jtfs(n = 2048, J = 6, J_fr = 3, Q = 8),
  "jtfs joint list" = function() scattering_jtfs(n = 2048, J = 6, J_fr = 3, Q = 8, format = "joint", out_type = "list"),
  "jtfs joint array" = function() scattering_jtfs(n = 2048, J = 6, J_fr = 3, Q = 8, format = "joint", F = 4),
  "jtfs F global" = function() scattering_jtfs(n = 2048, J = 6, J_fr = 3, Q = 8, F = "global"),
  "jtfs F 0" = function() scattering_jtfs(n = 2048, J = 6, J_fr = 3, Q = 8, F = 0, out_type = "list"),
  "jtfs T global" = function() scattering_jtfs(n = 2048, J = 6, J_fr = 3, Q = 8, T = "global"),
  "jtfs video" = function() scattering_jtfs(n = 900, J = 7, J_fr = 3, Q = c(8, 1), T_sec = 6, sr = 30)
)

test_that("the Rust and R engines agree on every operator", {
  for (name in names(operators)) {
    op <- operators[[name]]()
    x <- signal(op$n)
    a <- with_engine("r", scat_transform(op, x))
    b <- with_engine("rust", scat_transform(op, x))
    expect_identical(b$meta, a$meta, label = name)
    expect_identical(b$spec$engine, "rust")
    expect_identical(a$spec$engine, "r")
    expect_equal(b$coef, a$coef, tolerance = 1e-10, label = name)
  }
})

test_that("joint coefficients renormalise in every storage, and lose amplitude", {
  for (name in grep("^jtfs", names(operators), value = TRUE)) {
    op <- operators[[name]]()
    if (!identical(op$average, "local")) next
    x <- signal(op$n)
    a <- with_engine("rust", scat_renorm(scat_transform(op, x), eps = 1e-300))
    b <- with_engine("rust", scat_renorm(scat_transform(op, 40 * x), eps = 1e-300))
    o2 <- which(a$meta$order == 2L)
    pick <- function(s) if (is.list(s$coef)) unlist(s$coef[o2]) else as.numeric(asplit(s$coef, 1L)[o2] |> unlist())
    # Joint arrays keep padded frequency rows beyond the real bands, where
    # path and denominator are both tiny and the ratio loses digits; within
    # the real bands agreement is about 1e-9.
    expect_equal(pick(b), pick(a), tolerance = 1e-6, label = name)
    expect_true(all(is.finite(pick(a))), label = name)
  }
})

test_that("results record the engine and its numerics", {
  op <- scattering_1d(n = 1024, J = 5)
  s <- with_engine("rust", scat_transform(op, signal(1024)))
  expect_identical(s$spec$numerics, wavscatengine::engine_version()$numerics)
  expect_true(is.na(with_engine("r", scat_transform(op, signal(1024)))$spec$numerics))
})

test_that("the engine reproduces the reference bits on this machine", {
  expect_length(wavscatengine::engine_verify(), 0L)
})

test_that("engine choice is validated", {
  expect_error(with_engine("fortran", scat_engine()), "must be")
  expect_identical(with_engine("r", scat_engine()), "r")
  expect_identical(with_engine("auto", scat_engine()), "rust")
})
