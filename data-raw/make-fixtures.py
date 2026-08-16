"""Generate Kymatio reference values for the wavscat parity tests.

Development-time only. Kymatio is never a dependency of the package; this
script writes raw float64 arrays plus a manifest, which data-raw/make-fixtures.R
then converts into tests/testthat/fixtures/*.rds.

Requires Kymatio at or after the `anden_generator` refactor (0.4.0.dev0 or
later), which also carries the fix to issue #984 in the second-order energy
scaling. Install with:

    pip install --no-deps git+https://github.com/kymatio/kymatio.git@main

Usage:

    python data-raw/make-fixtures.py <output-directory>
"""
import os
import sys
import math

import numpy as np

from kymatio.scattering1d.filter_bank import (
    adaptive_choice_P, morlet_1d, gauss_1d, compute_sigma_psi,
    compute_temporal_support, get_max_dyadic_subsampling, compute_xi_max,
    anden_generator, scattering_filter_factory)
from kymatio.scattering1d.frontend.numpy_frontend import ScatteringNumPy1D

R_PSI = math.sqrt(0.5)
SIGMA0 = 0.1
ALPHA = 5.0

_manifest = []
_outdir = None


def put(name, arr):
    a = np.asarray(arr, dtype=np.float64).ravel()
    a.tofile(os.path.join(_outdir, name + ".bin"))
    _manifest.append((name, len(a)))


# --------------------------------------------------------------------------
# Filter bank internals
# --------------------------------------------------------------------------

FB_SIGMAS = [0.001, 0.01, 0.0208, 0.05, 0.1, 0.2, 0.5]
FB_QS = [1, 2, 4, 8, 16]
FB_PAIRS = [(0.4, 0.02), (0.1, 0.01), (0.01, 0.001), (0.35, 0.15),
            (0.001, 0.0001)]
FB_MORLET = [(256, 0.4, 0.0208), (256, 0.05, 0.003), (512, 0.35, 0.15),
             (1024, 0.01, 0.0005), (128, 0.3, 0.25)]
FB_GAUSS = [(256, 0.0125), (512, 0.05), (1024, 0.000390625), (128, 0.3)]
FB_TSUPP = [(1024, 0.0125), (2048, 0.00390625), (512, 0.05)]
FB_GEN = [(6, 1), (6, 8), (8, 8), (10, 2), (4, 16), (12, 1)]
FB_FACTORY = [(2048, 6, (8, 1), 64), (1024, 5, (4, 2), 32),
              (4096, 8, (8, 2), 256), (512, 4, (1, 1), 16)]


def dump_filter_bank():
    put("scalar_sigmas", FB_SIGMAS)
    put("scalar_P", [adaptive_choice_P(s) for s in FB_SIGMAS])
    put("scalar_Q", FB_QS)
    put("scalar_xi_max", [compute_xi_max(q) for q in FB_QS])
    put("scalar_sigma_psi", [compute_sigma_psi(0.4, q, R_PSI) for q in FB_QS])
    put("scalar_pair_xi", [p[0] for p in FB_PAIRS])
    put("scalar_pair_sigma", [p[1] for p in FB_PAIRS])
    put("scalar_maxsub",
        [get_max_dyadic_subsampling(x, s, ALPHA) for x, s in FB_PAIRS])

    put("morlet_meta", [v for c in FB_MORLET for v in c])
    for i, (N, xi, sig) in enumerate(FB_MORLET):
        put("morlet_%d" % i, morlet_1d(N, xi, sig))

    put("gauss_meta", [v for c in FB_GAUSS for v in c])
    for i, (N, sig) in enumerate(FB_GAUSS):
        put("gauss_%d" % i, gauss_1d(N, sig))

    put("tsupp_meta", [v for c in FB_TSUPP for v in c])
    put("tsupp", [compute_temporal_support(gauss_1d(N, s).reshape(1, -1))
                  for N, s in FB_TSUPP])

    put("gen_meta", [v for c in FB_GEN for v in c])
    for i, (J, Q) in enumerate(FB_GEN):
        spec = list(anden_generator(J, Q, sigma0=SIGMA0, r_psi=R_PSI))
        put("gen_xi_%d" % i, [s[0] for s in spec])
        put("gen_sigma_%d" % i, [s[1] for s in spec])

    fb_kwargs = {"alpha": ALPHA, "r_psi": R_PSI, "sigma0": SIGMA0}
    put("fac_meta",
        [v for N, J, Q, T in FB_FACTORY for v in (N, J, Q[0], Q[1], T)])
    for i, (N, J, Q, T) in enumerate(FB_FACTORY):
        phi, psi1, psi2 = scattering_filter_factory(
            N, J, Q, T, (anden_generator, fb_kwargs), np.sum)
        put("fac%d_phi_j" % i, [phi["j"]])
        put("fac%d_phi_sigma" % i, [phi["sigma"]])
        put("fac%d_phi_nlev" % i, [len(phi["levels"])])
        for lev, v in enumerate(phi["levels"]):
            put("fac%d_phi_lev%d" % (i, lev), v)
        for order, bank in ((1, psi1), (2, psi2)):
            put("fac%d_p%d_n" % (i, order), [len(bank)])
            put("fac%d_p%d_xi" % (i, order), [p["xi"] for p in bank])
            put("fac%d_p%d_sigma" % (i, order), [p["sigma"] for p in bank])
            put("fac%d_p%d_j" % (i, order), [p["j"] for p in bank])
            put("fac%d_p%d_nlev" % (i, order),
                [len(p["levels"]) for p in bank])
            for n, p in enumerate(bank):
                for lev, v in enumerate(p["levels"]):
                    put("fac%d_p%d_%d_lev%d" % (i, order, n, lev), v)


# --------------------------------------------------------------------------
# Full transforms
# --------------------------------------------------------------------------

# (N, J, Q, T, max_order, stride, out_type, signal)
CASES = [
    (2048, 6, 8,       None,     2, None, "array", "noise"),
    (2048, 6, 8,       None,     2, None, "array", "chirp"),
    (1024, 5, (4, 2),  32,       2, None, "array", "noise"),
    (4096, 8, (8, 2),  None,     2, None, "array", "am"),
    (2048, 6, 1,       None,     1, None, "array", "noise"),
    (2048, 6, 8,       None,     2, 16,   "array", "tone"),
    (2048, 6, 8,       "global", 2, None, "array", "noise"),
    (512,  4, 1,       None,     2, None, "array", "impulse"),
    (3000, 7, (12, 1), 100,      2, None, "array", "noise"),
    (2048, 6, 8,       0,        2, None, "list",  "noise"),
    (1500, 6, (8, 1),  30,       2, None, "array", "am"),
    (8192, 9, (8, 2),  None,     2, None, "array", "chirp"),
]


def signal(kind, N):
    t = np.arange(N) / N
    if kind == "noise":
        return np.random.RandomState(0).randn(N)
    if kind == "tone":
        return np.sin(2 * np.pi * 0.11 * np.arange(N))
    if kind == "chirp":
        return np.sin(2 * np.pi * (0.005 + 0.4 * t) * np.arange(N))
    if kind == "impulse":
        x = np.zeros(N)
        x[N // 3] = 1.0
        return x
    if kind == "am":
        carrier = np.sin(2 * np.pi * 0.15 * np.arange(N))
        envelope = 1 + 0.8 * np.sin(2 * np.pi * 0.004 * np.arange(N))
        return carrier * envelope
    raise ValueError(kind)


def dump_transforms():
    put("n_cases", [len(CASES)])
    for i, (N, J, Q, T, mo, stride, ot, kind) in enumerate(CASES):
        x = signal(kind, N)
        qq = Q if isinstance(Q, tuple) else (Q, 1)
        put("case%d_x" % i, x)
        put("case%d_cfg" % i,
            [N, J, qq[0], qq[1], mo, -1 if stride is None else stride])
        # T encoding: -1 for None, -2 for "global", otherwise the value.
        put("case%d_T" % i,
            [-1 if T is None else (-2 if T == "global" else T)])
        put("case%d_out_type" % i, [0 if ot == "array" else 1])

        s = ScatteringNumPy1D(J=J, shape=N, Q=Q, T=T, max_order=mo,
                              stride=stride, out_type=ot)
        out = s(x)
        meta = s.meta()
        put("case%d_meta_order" % i, meta["order"])
        put("case%d_meta_n" % i, np.nan_to_num(meta["n"], nan=-1.0))
        put("case%d_meta_j" % i, np.nan_to_num(meta["j"], nan=-1.0))
        put("case%d_meta_xi" % i, np.nan_to_num(meta["xi"], nan=-1.0))

        if ot == "array":
            put("case%d_dim" % i, list(out.shape))
            put("case%d_out" % i, out)
        else:
            put("case%d_npaths" % i, [len(out)])
            put("case%d_lens" % i, [len(p["coef"]) for p in out])
            for k, p in enumerate(out):
                put("case%d_p%d" % (i, k), p["coef"])


def main():
    global _outdir
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    _outdir = sys.argv[1]
    os.makedirs(_outdir, exist_ok=True)
    dump_filter_bank()
    dump_transforms()
    with open(os.path.join(_outdir, "manifest.tsv"), "w") as fh:
        for name, n in _manifest:
            fh.write("%s\t%d\n" % (name, n))
    print("wrote %d arrays to %s" % (len(_manifest), _outdir))


if __name__ == "__main__":
    main()
