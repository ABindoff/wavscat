"""Generate Kymatio reference values for the wavscat JTFS parity tests.

Development-time only; see make-fixtures.py for installation notes.

Coefficients are dumped from Kymatio's `out_type="dict"`, keyed by path, rather
than from its stacked array. That is deliberate. As of 0.4.0.dev0, Kymatio's
`TimeFrequencyScattering` stacks `format="time"` coefficients in generator order
(all bands of one frequential filter, then the next) while its `meta()` reports
them in sorted order (all filters of one band, then the next), so the array and
its metadata disagree. `Scattering1D` sorts before stacking; the joint frontend
does not. The dict output is keyed by path and therefore unambiguous, so it is
what we compare against.

Usage:

    python data-raw/make-fixtures-jtfs.py <output-directory>
"""
import os
import sys

import numpy as np

from kymatio.scattering1d.frontend.numpy_frontend import (
    TimeFrequencyScatteringNumPy)

_manifest = []
_outdir = None


def put(name, arr):
    a = np.asarray(arr, dtype=np.float64).ravel()
    a.tofile(os.path.join(_outdir, name + ".bin"))
    _manifest.append((name, len(a)))


def signal(kind, N):
    t = np.arange(N) / N
    if kind == "noise":
        return np.random.RandomState(0).randn(N)
    if kind == "chirp":
        return np.sin(2 * np.pi * (0.005 + 0.4 * t) * np.arange(N))
    if kind == "am":
        carrier = np.sin(2 * np.pi * 0.15 * np.arange(N))
        envelope = 1 + 0.8 * np.sin(2 * np.pi * 0.004 * np.arange(N))
        return carrier * envelope
    raise ValueError(kind)


# (N, J, J_fr, Q, Q_fr, T, F, format, signal)
CASES = [
    (4096, 6, 3, 8,      1, None,     None,     "time",  "noise"),
    (4096, 6, 3, 8,      1, None,     None,     "joint", "chirp"),
    (2048, 5, 2, 4,      1, None,     None,     "time",  "am"),
    (4096, 6, 3, (8, 2), 1, None,     None,     "time",  "chirp"),
    (4096, 6, 3, 8,      2, None,     None,     "time",  "noise"),
    (4096, 6, 3, 8,      1, "global", None,     "time",  "noise"),
    (4096, 6, 3, 8,      1, None,     "global", "time",  "chirp"),
    (4096, 6, 3, 8,      1, None,     4,        "joint", "am"),
    (8192, 7, 4, 8,      1, None,     None,     "time",  "chirp"),
    (4096, 6, 3, 8,      1, "global", "global", "joint", "noise"),
    (4096, 6, 3, 8,      1, None,     0,        "time",  "am"),
    (4096, 6, 3, 8,      1, 256,      8,        "time",  "chirp"),
]


def main():
    global _outdir
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    _outdir = sys.argv[1]
    os.makedirs(_outdir, exist_ok=True)

    put("n_cases", [len(CASES)])
    for i, (N, J, J_fr, Q, Q_fr, T, F, fmt, kind) in enumerate(CASES):
        x = signal(kind, N)
        qq = Q if isinstance(Q, tuple) else (Q, 1)
        put("case%d_x" % i, x)
        put("case%d_cfg" % i, [N, J, J_fr, qq[0], qq[1], Q_fr,
                               0 if fmt == "time" else 1])
        # T and F encoding: -1 for None, -2 for "global", else the value.
        put("case%d_T" % i,
            [-1 if T is None else (-2 if T == "global" else T)])
        put("case%d_F" % i,
            [-1 if F is None else (-2 if F == "global" else F)])

        s = TimeFrequencyScatteringNumPy(
            J=J, J_fr=J_fr, Q=Q, Q_fr=Q_fr, shape=N, T=T, F=F,
            format=fmt, out_type="dict")
        out = s(x)
        keys = list(out.keys())

        put("case%d_npaths" % i, [len(keys)])
        # Keys are ragged tuples, so pad to three columns with -1.
        put("case%d_keys" % i,
            [list(k) + [-1] * (3 - len(k)) for k in keys])
        put("case%d_keylen" % i, [len(k) for k in keys])
        put("case%d_ndim" % i, [np.asarray(out[k]).ndim for k in keys])
        put("case%d_shape" % i,
            [list(np.asarray(out[k]).shape) + [-1] *
             (2 - np.asarray(out[k]).ndim) for k in keys])
        for p, k in enumerate(keys):
            put("case%d_c%d" % (i, p), out[k])

    with open(os.path.join(_outdir, "manifest.tsv"), "w") as fh:
        for name, n in _manifest:
            fh.write("%s\t%d\n" % (name, n))
    print("wrote %d arrays for %d cases" % (len(_manifest), len(CASES)))


if __name__ == "__main__":
    main()
