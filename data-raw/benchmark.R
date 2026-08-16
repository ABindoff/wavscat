## Timings for the two target workloads, plus a few reference points.
##
##   Rscript data-raw/benchmark.R
##
## Operator construction is reported separately from application, because an
## operator is built once and applied to every recording in a study.

suppressMessages(devtools::load_all(".", quiet = TRUE))

time_it <- function(expr, reps = 5) {
  expr <- substitute(expr)
  env <- parent.frame()
  eval(expr, env)  # warm up
  t <- vapply(seq_len(reps), function(i) {
    system.time(eval(expr, env))[["elapsed"]]
  }, 0)
  stats::median(t)
}

row <- function(label, build, apply_, n_paths, n_time) {
  cat(sprintf("%-44s  build %7.3f s   apply %7.3f s   %5d paths x %s\n",
              label, build, apply_, n_paths, n_time))
}

cat("R", R.version.string, "\n\n")
cat("=== target workloads ===\n")

## Finger tapping: 60 s of webcam landmarks at 30 fps -----------------------
set.seed(1)
tap <- as.numeric(stats::filter(stats::rnorm(1800), 0.8, "recursive"))
b <- time_it(scat_preset("kinematic", n = 1800, sr = 30))
op <- scat_preset("kinematic", n = 1800, sr = 30)
a <- time_it(scat_transform(op, tap), reps = 20)
row("tapping, 1800 samples @ 30 Hz", b, a, nrow(op$paths),
    dim(scat_transform(op, tap)$coef)[2])

## Ten landmark channels at once -------------------------------------------
taps <- matrix(stats::rnorm(1800 * 10), nrow = 1800)
a <- time_it(scat_transform(op, taps), reps = 5)
row("  same, 10 channels", 0, a, nrow(op$paths), "-")

## DDK audio: 10 s at 44.1 kHz ---------------------------------------------
n_audio <- 441000L
audio <- as.numeric(stats::filter(stats::rnorm(n_audio), 0.5, "recursive"))
b <- time_it(scat_preset("audio", n = n_audio, sr = 44100), reps = 3)
op_a <- scat_preset("audio", n = n_audio, sr = 44100)
a <- time_it(scat_transform(op_a, audio), reps = 3)
row("DDK audio, 10 s @ 44.1 kHz", b, a, nrow(op_a$paths),
    dim(scat_transform(op_a, audio)$coef)[2])

## DDK audio at 16 kHz, the usual speech rate -------------------------------
n16 <- 160000L
a16 <- as.numeric(stats::filter(stats::rnorm(n16), 0.5, "recursive"))
b <- time_it(scat_preset("speech", n = n16, sr = 16000), reps = 3)
op_s <- scat_preset("speech", n = n16, sr = 16000)
a <- time_it(scat_transform(op_s, a16), reps = 3)
row("DDK audio, 10 s @ 16 kHz", b, a, nrow(op_s$paths),
    dim(scat_transform(op_s, a16)$coef)[2])

cat("\n=== joint time-frequency scattering ===\n")
op_j <- scattering_jtfs(n = 32768, J = 10, J_fr = 3, Q = 8, sr = 16000)
xj <- as.numeric(stats::filter(stats::rnorm(32768), 0.5, "recursive"))
b <- time_it(scattering_jtfs(n = 32768, J = 10, J_fr = 3, Q = 8, sr = 16000), reps = 3)
a <- time_it(scat_transform(op_j, xj), reps = 3)
row("JTFS, 2 s @ 16 kHz (32768 samples)", b, a, nrow(op_j$paths), "-")

cat("\n=== scaling with signal length (J = 8, Q = 8) ===\n")
for (n in c(2^12, 2^14, 2^16, 2^18, 2^20)) {
  op_n <- scattering_1d(n = n, J = 8, Q = 8)
  xn <- stats::rnorm(n)
  b <- time_it(scattering_1d(n = n, J = 8, Q = 8), reps = 3)
  a <- time_it(scat_transform(op_n, xn), reps = 3)
  row(sprintf("  n = 2^%d", log2(n)), b, a, nrow(op_n$paths), "-")
}

cat("\n=== a realistic batch: 200 tapping trials ===\n")
batch <- matrix(stats::rnorm(1800 * 200), nrow = 1800)
el <- system.time({
  S <- scat_transform(op, batch)
  f <- scat_features(scat_log(scat_renorm(S), eps = scat_eps(S)))
})[["elapsed"]]
cat(sprintf("  transform + renorm + log + features : %.2f s  (%.1f ms per trial)\n",
            el, 1000 * el / 200))
cat(sprintf("  feature table: %d rows x %d columns\n", nrow(f), ncol(f)))
