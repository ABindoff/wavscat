# wavscat

Wavelet scattering transforms for time series and audio, in base R.

Wavelet scattering builds a signal representation that is stable to deformation
and locally invariant to translation, while keeping the high-frequency detail
that a spectrogram averages away. Nothing is learned from data, so it behaves
well on the small, noisy, variable-timing recordings common in behavioural work.

`wavscat` implements 1D time scattering of order 0, 1 and 2, and joint
time-frequency scattering. Everything is computed in the Fourier domain using
base R, so there is nothing to compile.

## Installation

```r
# install.packages("pak")
pak::pak("ABindoff/wavscat")
```

## Usage

```r
library(wavscat)

# Build an operator once, apply it to as many signals as you like.
op <- scattering_1d(n = 4096, J = 8, Q = 8, sr = 8000)
S  <- scat_transform(op, x)

# Coefficients, plus a table saying what each one means.
S$coef            # [path, time, channel]
scat_meta(S)      # order, centre frequencies in Hz, dyadic scales

# Model-ready features.
features <- scat_features(scat_log(scat_renorm(S), eps = scat_eps(S)))
```

Presets carry sensible parameters for two common kinds of recording:

```r
scat_preset("kinematic", n = 1800, sr = 30)      # webcam-rate movement
scat_preset("speech",    n = 32000, sr = 16000)  # speech and DDK audio
```

Joint time-frequency scattering adds a filter bank along the log-frequency axis,
which makes the representation sensitive to spectral energy drifting up or down
in frequency. A rising formant transition and a falling one have identical time
scattering coefficients and opposite *spins* here:

```r
op <- scattering_jtfs(n = 32000, J = 10, J_fr = 3, Q = 8, sr = 16000)
```

Inside tidymodels, `step_scattering()` learns the log-compression scale from the
training half of each resample, which doing the transform by hand before
splitting gets wrong.

## Correctness

The filter bank and both cascades are checked against
[Kymatio](https://github.com/kymatio/kymatio), the reference implementation, to
a relative error of about `1e-15`: every filter at every resolution, and 24
transform configurations spanning both formats, all averaging modes, custom
strides, and non-power-of-two signal lengths. Kymatio is a development-time
dependency only; the fixtures are committed and the test suite needs no Python.

A separate set of tests checks properties that a faithful transcription of a
misread formula would still satisfy: the Littlewood-Paley sum staying at its
frame bound of one, wavelets having exactly zero mean, homogeneity,
non-expansiveness, energy decaying with order, and translation sensitivity
growing linearly in the shift over `T`.

Two Kymatio behaviours are deliberately not reproduced, both documented in the
test suite. Its `format = "time"` joint output is stacked in generator order
while its `meta()` is sorted, so the two disagree; `wavscat` sorts both.

## Performance

Single recording, ordinary laptop, base R:

| workload | time |
|---|---|
| tapping trace, 1800 samples at 30 Hz | ~20 ms |
| DDK audio, 10 s at 16 kHz | ~1.7 s |
| DDK audio, 10 s at 44.1 kHz | ~4.2 s |
| joint scattering, 2 s at 16 kHz | ~5.5 s |

Building an operator costs a few seconds at audio lengths, but happens once and
is reused across every recording of the same length. For audio-scale work,
parallelise across recordings; each is independent.

## Optional Rust engine

wavscat needs no compiler. If you also install the companion package
[wavscatengine](https://github.com/ABindoff/wavscat-core/tree/main/r/wavscatengine),
wavscat uses it automatically: the same Rust code that runs in Python and in
the browser, so coefficients and features are bit-for-bit identical across
languages and platforms, and joint scattering is about four times faster
(the 2 s, 16 kHz example above takes about 1.3 s instead of 5.7 s).

```r
# Needs Rust (https://rustup.rs) and, on Windows, Rtools and
# `rustup target add x86_64-pc-windows-gnu`.
pak::pak("ABindoff/wavscat-core/r/wavscatengine")
scat_engine()   # "rust" once installed
```

The two engines agree to about 1e-14 in any case. `options(wavscat.engine =
"r")` keeps the base R engine; each result records its engine in
`x$spec$engine`. Renormalising joint time-frequency coefficients with
`scat_renorm()` needs the Rust engine.

## References

Andén, J. and Mallat, S. (2014). Deep scattering spectrum. *IEEE Transactions on
Signal Processing*, 62(16), 4114–4128.

Andén, J., Lostanlen, V. and Mallat, S. (2019). Joint time-frequency scattering.
*IEEE Transactions on Signal Processing*, 67(14), 3704–3718.

Andreux, M. et al. (2020). Kymatio: scattering transforms in Python. *Journal of
Machine Learning Research*, 21(60), 1–6.
