make_S <- function(n = 2048, ...) {
  set.seed(5)
  scattering(cbind(a = stats::rnorm(n), b = stats::rnorm(n)), J = 6, Q = 8, ...)
}


test_that("log compression is monotone and preserves shape", {
  S <- make_S()
  L <- scat_log(S, eps = 1e-3)
  expect_equal(dim(L$coef), dim(S$coef))
  expect_true(all(is.finite(L$coef)))
  # Monotone, so the ranking of coefficients within a channel is unchanged.
  expect_equal(order(L$coef[, 1, 1]), order(S$coef[, 1, 1]))
  expect_equal(L$transforms, "log")

  # Orders one and up are non-negative, and there it is plain log1p.
  up <- S$meta$order >= 1L
  expect_true(all(L$coef[up, , ] >= 0))
  expect_equal(L$coef[up, , ], log1p(S$coef[up, , ] / 1e-3))
})


test_that("log compression handles the signed order-zero coefficient", {
  S <- make_S()
  # The local mean genuinely goes negative, which is why plain log1p is unsafe.
  s0 <- S$coef[S$meta$order == 0L, , ]
  expect_true(any(s0 < 0))

  L <- scat_log(S, eps = 1e-3)
  l0 <- L$coef[L$meta$order == 0L, , ]
  expect_true(all(is.finite(l0)))
  expect_equal(sign(l0), sign(s0))
  expect_equal(l0, sign(s0) * log1p(abs(s0) / 1e-3))
})


test_that("log compression validates eps and warns on reapplication", {
  S <- make_S()
  expect_error(scat_log(S, eps = 0), "positive number")
  expect_error(scat_log(S, eps = c(1, 2)), "positive number")
  expect_warning(scat_log(scat_log(S)), "already been log-compressed")
  expect_error(scat_log(1:10), "wavscat_coefs")
})


test_that("scat_eps sits inside the range of the coefficients", {
  S <- make_S()
  e <- scat_eps(S)
  vals <- S$coef[S$meta$order >= 1, , ]
  expect_gt(e, 0)
  expect_gt(e, min(vals))
  expect_lt(e, stats::median(vals))
  # A list of objects is pooled.
  expect_gt(scat_eps(list(S, S)), 0)
})


test_that("renormalisation divides each order-2 path by its parent", {
  S <- make_S()
  R <- scat_renorm(S, eps = 0)
  meta <- S$meta
  o2 <- which(meta$order == 2L)
  o1 <- which(meta$order == 1L)

  for (k in c(1L, length(o2) %/% 2L, length(o2))) {
    row <- o2[k]
    parent <- o1[match(meta$n1[row], meta$n1[o1])]
    expect_equal(R$coef[row, , ], S$coef[row, , ] / S$coef[parent, , ])
  }
  # Orders 0 and 1 are untouched.
  expect_equal(R$coef[meta$order < 2L, , ], S$coef[meta$order < 2L, , ])
})


test_that("renormalisation must come before log compression", {
  S <- make_S()
  expect_error(scat_renorm(scat_log(S)), "Renormalise before log")
  expect_warning(scat_renorm(scat_renorm(S)), "already been renormalised")
})


test_that("renormalisation removes the effect of overall amplitude on order 2", {
  # The point of the step: a louder recording of the same modulation should give
  # the same renormalised order-2 coefficients.
  set.seed(9)
  op <- scattering_1d(n = 4096, J = 8, Q = 8)
  carrier <- sin(2 * pi * 0.15 * seq_len(4096))
  x <- carrier * (1 + 0.8 * sin(2 * pi * 0.004 * seq_len(4096)))
  o2 <- scat_meta(op)$order == 2L

  # Exactly invariant only with eps = 0; any positive eps deliberately breaks it
  # for bands that carry no energy.
  a <- scat_renorm(scat_transform(op, x), eps = 0)$coef[o2, , 1]
  b <- scat_renorm(scat_transform(op, 25 * x), eps = 0)$coef[o2, , 1]
  expect_lt(max(abs(a - b)) / max(abs(a)), 1e-9)
})


test_that("features come out one row per channel with decodable names", {
  S <- make_S()
  f <- scat_features(S)
  expect_s3_class(f, "tbl_df")
  expect_equal(nrow(f), 2L)
  expect_equal(f$channel, c("a", "b"))
  expect_equal(ncol(f), nrow(S$meta) + 1L)
  expect_equal(names(f)[-1], S$meta$path)
  expect_true(all(vapply(f[-1], is.numeric, TRUE)))
  expect_false(anyNA(f))

  # The values really are the per-path time means.
  expect_equal(as.numeric(f[1, -1]), unname(apply(S$coef[, , 1], 1, mean)))
})


test_that("feature summaries and order selection behave", {
  S <- make_S()
  expect_equal(
    as.numeric(scat_features(S, summary = "max")[1, -1]),
    unname(apply(S$coef[, , 1], 1, max))
  )
  expect_equal(
    as.numeric(scat_features(S, summary = "sd")[1, -1]),
    unname(apply(S$coef[, , 1], 1, stats::sd))
  )

  f1 <- scat_features(S, orders = 1)
  expect_equal(ncol(f1) - 1L, sum(S$meta$order == 1L))
  expect_true(all(grepl("^S1_", names(f1)[-1])))
  expect_error(scat_features(S, orders = 5), "No paths of order")

  f_pre <- scat_features(S, prefix = "tap.")
  expect_true(all(grepl("^tap\\.", names(f_pre)[-1])))
})


test_that("summary = 'none' keeps the time axis in the right order", {
  S <- make_S()
  n_time <- dim(S$coef)[2L]
  f <- scat_features(S, summary = "none", orders = 1)
  n_paths <- sum(S$meta$order == 1L)

  expect_equal(ncol(f) - 1L, n_paths * n_time)
  # Time varies fastest, so each path's course is contiguous.
  expect_equal(names(f)[2:(n_time + 1)],
               paste0(S$meta$path[S$meta$order == 1L][1], "_t", seq_len(n_time)))
  expect_equal(as.numeric(f[1, 2:(n_time + 1)]),
               unname(S$coef[which(S$meta$order == 1L)[1], , 1]))
})


test_that("feature names are syntactically valid", {
  f <- scat_features(make_S())
  expect_equal(names(f), make.names(names(f), unique = TRUE))
})


test_that("the recommended pipeline runs end to end", {
  S <- make_S()
  out <- scat_features(scat_log(scat_renorm(S), eps = scat_eps(S)))
  expect_equal(nrow(out), 2L)
  expect_true(all(is.finite(as.matrix(out[-1]))))
})
