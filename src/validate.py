"""Validate our reading of the paper's conventions against Supplementary Data 1.

Part 1 (conventions, all 344 recordings):
  1. Total 1-50 Hz == 0.5 * sum of all 98 bins.
  2. Delta 2-4 Hz: which integration rule reproduces the reported Fig 1b delta,
     and whether the same rule applied to Fig 1c reproduces Fig 1d periodic delta.
  3. log10(Fig 1a) - log10(Fig 1e) vs Fig 1c.

Band rule found (Nov 1): band power = (hi - lo) * mean(psd at bins lo <= f <= hi),
i.e. 0.4 * sum of the five bins 2.0..4.0 Hz for delta.

Writes data/validate_delta.csv (per-recording residuals, local only: CC BY-NC-ND)
and tables/delta_rules.csv (aggregate summary, published).

Usage: python src/validate.py [path/to/supplementary-data-1.xlsx]
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

import load_supplement as ls
from fit_spectra import DELTA, band_power

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_XLSX = ROOT / "refs" / "supplementary-data-1.xlsx"
DATA = ROOT / "data"  # per-recording derived values: local only (CC BY-NC-ND), gitignored
TABLES = ROOT / "tables"  # aggregate summaries: published
DF = 0.5  # Hz, bin width


def delta_rules(freqs, psd):
    """Candidate 2-4 Hz integrals, each shape (n_recordings,)."""
    half_open = (freqs >= 2) & (freqs < 4)
    closed = (freqs >= 2) & (freqs <= 4)
    return {
        "sum[2,4)": psd[:, half_open].sum(axis=1),
        "sum[2,4]": psd[:, closed].sum(axis=1),
        "0.5*sum[2,4)": DF * psd[:, half_open].sum(axis=1),
        "0.5*sum[2,4]": DF * psd[:, closed].sum(axis=1),
        "trapezoid[2,4]": np.trapezoid(psd[:, closed], freqs[closed], axis=1),
        "2*mean[2,4]": band_power(freqs, psd, *DELTA),
    }


def check_total(psd, reported):
    ours = DF * psd.sum(axis=1)
    rel = ours / reported["total"].to_numpy() - 1
    bad = np.flatnonzero(~np.isclose(ours, reported["total"], rtol=1e-5))
    print("1. Total 1-50 Hz == 0.5 * sum(all 98 bins)  (= 49 Hz * mean)")
    print(f"   within rtol 1e-5: {rel.size - bad.size}/{rel.size}")
    for i in bad:
        print(f"   FAIL {reported.index[i]}: reported {reported['total'].iloc[i]:g}, "
              f"ours {ours[i]:g} ({rel[i]:+.2%})")
    return rel


def check_delta(freqs, psd, reported):
    target = reported["delta"].to_numpy()
    rules = delta_rules(freqs, psd)
    rel = {name: v / target - 1 for name, v in rules.items()}

    print("\n2a. Delta 2-4 Hz vs reported Fig 1b delta (relative error = ours / reported - 1)")
    print(f"   {'rule':16s}{'median':>11s}{'median|e|':>11s}{'p5':>10s}{'p95':>10s}{'max|e|':>10s}{'closest':>9s}")
    abs_err = np.vstack([np.abs(e) for e in rel.values()])
    closest = np.array(list(rules))[abs_err.argmin(axis=0)]
    for name, e in rel.items():
        q = np.percentile(e, [5, 95])
        print(f"   {name:16s}{np.median(e):11.2e}{np.median(np.abs(e)):11.2e}{q[0]:10.2e}{q[1]:10.2e}"
              f"{np.abs(e).max():10.2e}{(closest == name).sum():9d}")

    e = rel["trapezoid[2,4]"]
    print("   trapezoid[2,4] |relative error| distribution (rejected rule):")
    for lo, hi in [(0, 0.01), (0.01, 0.05), (0.05, 0.10), (0.10, np.inf)]:
        n = ((np.abs(e) >= lo) & (np.abs(e) < hi)).sum()
        print(f"     [{lo:.0%}, {hi:.0%}): {n}")

    summary = pd.DataFrame(
        {
            "median_rel_error": [np.median(e) for e in rel.values()],
            "median_abs_rel_error": [np.median(np.abs(e)) for e in rel.values()],
            "p5_rel_error": [np.percentile(e, 5) for e in rel.values()],
            "p95_rel_error": [np.percentile(e, 95) for e in rel.values()],
            "max_abs_rel_error": [np.abs(e).max() for e in rel.values()],
            "n_closest": [(closest == name).sum() for name in rules],
        },
        index=pd.Index(list(rules), name="rule"),
    )

    out = pd.DataFrame({"reported_delta": target, **rules}, index=reported.index)
    out["delta_rel_error"] = rel["2*mean[2,4]"]
    out["trapezoid_rel_error"] = e
    out["closest_rule"] = closest
    return out, summary


def check_periodic_delta(path, freqs, psd, reported):
    """Same band rule on Fig 1c (fitted peak model) and on the raw log residual."""
    _, fig1c, _ = ls.load_freq_sheet(path, "Fig 1c")
    _, fig1e, _ = ls.load_freq_sheet(path, "Fig 1e")
    target = reported["periodic_delta"].to_numpy()
    candidates = {
        "Fig 1c": band_power(freqs, fig1c, *DELTA),
        "log10(1a)-log10(1e)": band_power(freqs, np.log10(psd) - np.log10(fig1e), *DELTA),
    }
    print("\n2b. 2*mean[2,4] vs reported Fig 1d periodic delta (absolute difference, log10 units)")
    for name, v in candidates.items():
        d = np.abs(v - target)
        print(f"   {name:22s} median {np.median(d):.2e}  max {d.max():.2e}  "
              f"corr {np.corrcoef(v, target)[0, 1]:.4f}")
    return candidates["Fig 1c"] - target


def check_fig1c(path, freqs, psd, labels):
    _, fig1e, _ = ls.load_freq_sheet(path, "Fig 1e")
    _, fig1c, _ = ls.load_freq_sheet(path, "Fig 1c")
    if not (fig1e > 0).all():
        raise ValueError("Fig 1e: non-positive values")
    diff = np.abs((np.log10(psd) - np.log10(fig1e)) - fig1c)
    imax = diff.argmax(axis=1)
    f_at_max = freqs[imax]

    print("\n3. log10(Fig 1a) - log10(Fig 1e) vs Fig 1c (|difference|, log10 units)")
    print(f"   median over all bins and recordings = {np.median(diff):.4f}")
    k = diff.max(axis=1).argmax()
    print(f"   per-recording max: median = {np.median(diff.max(axis=1)):.4f}, "
          f"overall max = {diff.max():.4f} ({labels[k]} at {f_at_max[k]:g} Hz)")
    print(f"   frequency of per-recording max: median {np.median(f_at_max):g} Hz, "
          f"IQR {np.percentile(f_at_max, 25):g}-{np.percentile(f_at_max, 75):g} Hz")
    for lo, hi in [(1, 4), (4, 6), (6, 12), (12, 30), (30, 50)]:
        n = ((f_at_max >= lo) & (f_at_max < hi)).sum()
        print(f"     max in [{lo}, {hi}) Hz: {n}")
    return diff, f_at_max


def main(path=DEFAULT_XLSX):
    freqs, psd, labels = ls.load_spectra(path)
    reported = ls.load_reported(path)
    if list(reported.index) != labels:
        raise ValueError("Fig 1a and Fig 1b,d,f,g labels differ")

    total_rel = check_total(psd, reported)
    out, summary = check_delta(freqs, psd, reported)
    out["periodic_delta_abs_diff"] = check_periodic_delta(path, freqs, psd, reported)
    out["total_rel_error"] = total_rel
    diff, f_at_max = check_fig1c(path, freqs, psd, labels)
    out["fig1c_max_abs_diff"] = diff.max(axis=1)
    out["fig1c_max_abs_diff_hz"] = f_at_max

    DATA.mkdir(exist_ok=True)
    TABLES.mkdir(exist_ok=True)
    out.to_csv(DATA / "validate_delta.csv")
    summary.to_csv(TABLES / "delta_rules.csv", float_format="%.3g")
    print(f"\nwrote {DATA / 'validate_delta.csv'} (local only)")
    print(f"wrote {TABLES / 'delta_rules.csv'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_XLSX)
