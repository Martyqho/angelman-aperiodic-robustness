#!/usr/bin/env python3
"""Smoke test: refit the TD1 power spectrum from Supplementary Data 1 and
compare against the aperiodic parameters reported in the paper
(Communications Medicine, s43856-026-01741-0).

Reads:
  * 'Fig 1a'        -> Frequency + TD1 power spectrum
  * 'Fig 1b,d,f,g'  -> reported TD1 aperiodic exponent (Fig 1f) and offset (Fig 1g)

Fits with specparam (formerly FOOOF) in 'knee' mode (the paper's settings), then in 'fixed' mode.

Expected: knee exponent ~1.93-1.94, offset ~2.31, R^2 ~0.998; fixed exponent ~1.53.

Usage:
  pip install "specparam==2.0.0rc2" matplotlib openpyxl numpy
  python smoke_test.py [path/to/supplementary-data-1.xlsx]
"""

import sys

import numpy as np
import openpyxl
import specparam
from specparam import SpectralModel

DEFAULT_XLSX = "supplementary-data-1.xlsx"
SUBJECT = "TD1"

# Paper's settings (Methods, "EEG power spectral analysis"; specparam v2.0.0-rc2):
# "Fits were performed from 1 to 49.5 Hz using the 'knee' parameter ...
#  peak_width_limits = 0.5-12 Hz and peak_threshold = 2.0; ... max_n_peaks = inf".
# min_peak_height is not stated, so the library default (0) is used.
# Pinned to specparam 2.0.0rc2 to match the paper (later 2.0 RCs change the API).
FREQ_RANGE = [1, 49.5]
SPECPARAM_SETTINGS = dict(
    peak_width_limits=(0.5, 12),
    max_n_peaks=float("inf"),
    min_peak_height=0.0,
    peak_threshold=2.0,
)

# Expected refit values: name -> (aperiodic mode, result key, min, max), inclusive.
TARGETS = {
    "knee_exponent": ("knee", "exponent", 1.93, 1.94),
    "knee_offset": ("knee", "offset", 2.30, 2.32),
    "knee_r2": ("knee", "r2", 0.997, 1.0),
    "fixed_exponent": ("fixed", "exponent", 1.52, 1.54),
}


def read_spectrum(wb, subject):
    """Return (freqs, power) for `subject` from the 'Fig 1a' sheet."""
    ws = wb["Fig 1a"]
    rows = list(ws.iter_rows(values_only=True))
    # Row 1 holds group labels (TD / AS); row 2 holds 'Frequency', 'TD1', ...
    header_idx = next(i for i, r in enumerate(rows) if r and r[0] == "Frequency")
    header = rows[header_idx]
    col = header.index(subject)
    freqs, power = [], []
    for r in rows[header_idx + 1:]:
        if r[0] is None or r[col] is None:
            continue
        freqs.append(float(r[0]))
        power.append(float(r[col]))
    return np.asarray(freqs), np.asarray(power)


def read_reported(wb, subject):
    """Return (exponent, offset) for `subject` from the 'Fig 1b,d,f,g' sheet."""
    ws = wb["Fig 1b,d,f,g"]
    rows = list(ws.iter_rows(values_only=True))
    header = next(r for r in rows if r and subject in r)
    col = header.index(subject)
    exponent = offset = None
    for r in rows:
        label = str(r[1] or "").lower()
        if "aperiodic" in label and "exponent" in label:
            exponent = float(r[col])
        elif "aperiodic" in label and "offset" in label:
            offset = float(r[col])
    if exponent is None or offset is None:
        raise ValueError(f"Could not find reported aperiodic exponent/offset for {subject}")
    return exponent, offset


def fit(freqs, power, mode):
    fm = SpectralModel(aperiodic_mode=mode, verbose=False, **SPECPARAM_SETTINGS)
    fm.fit(freqs, power, FREQ_RANGE)
    ap = fm.aperiodic_params_
    # knee: [offset, knee, exponent]; fixed: [offset, exponent]
    result = {"offset": ap[0], "exponent": ap[-1], "r2": fm.r_squared_,
              "error": fm.error_, "n_peaks": len(fm.peak_params_)}
    if mode == "knee":
        # specparam fits the knee parameter k; the paper reports the knee
        # frequency F_knee = k^(1/exponent) in Hz.
        result["knee_k"] = ap[1]
        result["knee_freq"] = ap[1] ** (1 / ap[-1])
    return result


def check_targets(results):
    """Return [(name, value, lo, hi, ok)] for TARGETS, given {mode: fit result}."""
    rows = []
    for name, (mode, key, lo, hi) in TARGETS.items():
        value = results[mode][key]
        rows.append((name, value, lo, hi, lo <= value <= hi))
    return rows


def report(mode, rep_exp, rep_off, res):
    print(f"\n=== {SUBJECT} | aperiodic_mode='{mode}' | freq_range={FREQ_RANGE} ===")
    print(f"{'':12s}{'reported':>12s}{'refit':>12s}{'diff':>12s}")
    print(f"{'exponent':12s}{rep_exp:12.6f}{res['exponent']:12.6f}{res['exponent'] - rep_exp:+12.6f}")
    print(f"{'offset':12s}{rep_off:12.6f}{res['offset']:12.6f}{res['offset'] - rep_off:+12.6f}")
    if "knee_k" in res:
        print(f"{'knee k':12s}{'':>12s}{res['knee_k']:12.6f}   (model parameter, not a frequency)")
        print(f"{'knee freq':12s}{'':>12s}{res['knee_freq']:12.6f}   Hz = k^(1/exponent)")
    print(f"R^2 = {res['r2']:.4f}   error = {res['error']:.4f}   peaks = {res['n_peaks']}")


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_XLSX
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    freqs, power = read_spectrum(wb, SUBJECT)
    rep_exp, rep_off = read_reported(wb, SUBJECT)

    print(f"specparam {specparam.__version__}")
    print(f"Loaded {SUBJECT}: {len(freqs)} bins, {freqs[0]:g}-{freqs[-1]:g} Hz")
    print(f"specparam settings: {SPECPARAM_SETTINGS}")

    knee = fit(freqs, power, "knee")
    report("knee", rep_exp, rep_off, knee)

    fixed = fit(freqs, power, "fixed")
    report("fixed", rep_exp, rep_off, fixed)

    checks = check_targets({"knee": knee, "fixed": fixed})
    print("\n=== Targets ===")
    for name, value, lo, hi, ok in checks:
        print(f"[{'PASS' if ok else 'FAIL'}] {name:16s} {value:.4f} in [{lo}, {hi}]")
    sys.exit(0 if all(row[-1] for row in checks) else 1)


if __name__ == "__main__":
    main()
