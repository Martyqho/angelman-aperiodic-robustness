"""Secondary (a): agreement between the knee refit and the paper's reported values.

Reads data/features.csv (python src/fit_spectra.py) and the reported values from
Supplementary Data 1 ('Fig 1b,d,f,g'). For exponent, offset and periodic delta
(refit: band rule over the knee-mode fitted peak model), per group (TD, AS):
n, Pearson r, mean bias, SD, 95% limits of agreement (bias +- 1.96 SD), median
of refit - reported, median |refit - reported| and the fraction with |difference| > 0.15. 'all' includes every converged recording (the analysis set);
'no_knee_flag' excludes knee-flagged recordings (sensitivity).

Writes tables/validation.csv (aggregate) and tables/deviations.csv: every recorded
difference between this project's values and the paper's (AS159 Total 1-50 Hz, and
the knee-refit biases for exponent, offset and periodic delta, analysis set 'all').
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd

import load_supplement as ls
from config import load_config

ROOT = Path(__file__).resolve().parent.parent
FEATURES = {"exponent": ("exponent_knee", "exponent"),
            "offset": ("offset_knee", "offset"),
            "periodic_delta": ("periodic_delta_peak_knee", "periodic_delta")}


def agreement(refit, reported):
    d = refit - reported
    bias, sd = d.mean(), d.std(ddof=1)
    return {"n": len(d), "pearson_r": float(np.corrcoef(refit, reported)[0, 1]),
            "mean_refit_minus_reported": bias, "sd_diff": sd,
            "loa_lower": bias - 1.96 * sd, "loa_upper": bias + 1.96 * sd,
            "median_refit_minus_reported": float(np.median(d)),
            "median_abs_diff": float(np.median(np.abs(d))),
            "frac_abs_diff_gt_0.15": float((np.abs(d) > 0.15).mean()),
            "mean_reported": reported.mean(), "mean_refit": refit.mean()}


def main():
    logging.getLogger("load_supplement").setLevel(logging.ERROR)
    cfg = load_config()
    feats = pd.read_csv(ROOT / "data" / "features.csv").set_index("label")
    reported = ls.load_reported(ROOT / cfg["paper"]["supplement_file"])
    df = feats.join(reported[[r for _, r in FEATURES.values()]].add_suffix("_reported"))
    if len(df) != len(reported) or not df.index.equals(reported.index):
        raise RuntimeError("features.csv and reported values disagree on labels")
    df = df[df["success_knee"].astype(bool)]

    rows = []
    for subset, sub in [("all", df), ("no_knee_flag", df[~df["knee_flag"].astype(bool)])]:
        for g in ("TD", "AS"):
            s = sub[sub["group"] == g]
            for name, (ours, theirs) in FEATURES.items():
                col = f"{theirs}_reported"
                rows.append({"subset": subset, "group": g, "feature": name,
                             **agreement(s[ours].to_numpy(float), s[col].to_numpy(float))})
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "tables" / "validation.csv", index=False, float_format="%.5g")
    deviations(out, cfg).to_csv(ROOT / "tables" / "deviations.csv", index=False, float_format="%.5g")


def deviations(validation, cfg):
    """One row per recorded difference from the paper's reported values."""
    xlsx = ROOT / cfg["paper"]["supplement_file"]
    freqs, psd, labels = ls.load_spectra(xlsx)
    reported = ls.load_reported(xlsx)
    i = labels.index("AS159")
    ours_total = float(0.5 * psd[i].sum())
    rows = [{"item": "Total 1-50 Hz", "group": "AS159", "measure": "value (0.5 x sum of 98 bins)",
             "reported": reported.loc["AS159", "total"], "ours": ours_total,
             "difference": ours_total - reported.loc["AS159", "total"],
             "relative": ours_total / reported.loc["AS159", "total"] - 1,
             "source": "src/validate.py check 1; only recording of 344 off by > 1e-5"}]
    v = validation[validation["subset"] == "all"]
    for _, r in v.iterrows():
        rows.append({"item": f"{r['feature']}, knee refit", "group": r["group"],
                     "measure": "mean over recordings; difference = refit - reported",
                     "reported": r["mean_reported"], "ours": r["mean_refit"],
                     "difference": r["mean_refit_minus_reported"], "relative": float("nan"),
                     "source": f"tables/validation.csv (n {r['n']}, LoA {r['loa_lower']:.3f} to "
                               f"{r['loa_upper']:.3f}, median |d| {r['median_abs_diff']:.3f})"})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
