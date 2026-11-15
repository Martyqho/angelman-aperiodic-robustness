"""Follow-ups to the Fig 1e diagnostic (requested 2026-11-14, after the Nov 5 results; post hoc).

(i)  Subset test: per group, the periodic-delta gap (knee refit - reported) split by
     whether the refit's peak-model periodic delta is exactly 0 or > 0, with sizes. Added
     alongside (not requested): < 0.01 vs >= 0.01, because a fitted peak anywhere leaves a
     tiny Gaussian tail in 2-4 Hz, so "> 0" includes near-zero values.
(ii) Cross-tabs, 2 x 2 per group: Fig 1e status against the refit knee flag and
     against age < 3. Two status variables:
       unclassified  declared Fig 1e outcome 'unclassified' vs '(B)' (no (A)/(C) exist)
       flip          the Fig 1e single-curve residual changes sign between 1.0 and 1.5 Hz

Reads data/features.csv (src/fit_spectra.py) and data/fig1e_fits.csv
(src/fig1e_diagnostic.py). Writes tables/periodic_gap_subsets.csv and
tables/fig1e_crosstabs.csv (aggregate counts and summaries).
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import load_supplement as ls
from config import load_config
from fig1e_diagnostic import knee_curve

ROOT = Path(__file__).resolve().parent.parent
GROUPS = ("TD", "AS")


def gap_subsets(df):
    rows = []
    for g in GROUPS:
        sub = df[df["group"] == g]
        for name, m in (("refit_pd_eq_0", sub["pd_refit"] == 0), ("refit_pd_gt_0", sub["pd_refit"] > 0),
                        ("refit_pd_lt_0.01", sub["pd_refit"] < 0.01), ("refit_pd_ge_0.01", sub["pd_refit"] >= 0.01),
                        ("all", sub["pd_refit"].notna())):
            s = sub[m]
            d = s["pd_refit"] - s["pd_reported"]
            rows.append({"group": g, "subset": name, "n": len(s),
                         "mean_gap": d.mean(), "median_gap": d.median(),
                         "gap_q25": d.quantile(0.25), "gap_q75": d.quantile(0.75),
                         "mean_reported": s["pd_reported"].mean(), "mean_refit": s["pd_refit"].mean(),
                         "n_reported_eq_0": int((s["pd_reported"] == 0).sum()),
                         "share_of_group_gap": d.sum() / (sub["pd_refit"] - sub["pd_reported"]).sum()})
    return pd.DataFrame(rows)


def crosstabs(df):
    rows = []
    for g in GROUPS:
        sub = df[df["group"] == g]
        for status in ("unclassified", "flip"):
            for other in ("knee_flag", "age_lt_3"):
                t = pd.crosstab(sub[status], sub[other]).reindex(index=[True, False], columns=[True, False],
                                                                  fill_value=0)
                odds, p = stats.fisher_exact(t.to_numpy())
                rows.append({"group": g, "status": status, "versus": other,
                             "status_and_versus": int(t.loc[True, True]),
                             "status_not_versus": int(t.loc[True, False]),
                             "not_status_versus": int(t.loc[False, True]),
                             "not_status_not_versus": int(t.loc[False, False]),
                             "fisher_odds_ratio": odds, "fisher_p": p})
    return pd.DataFrame(rows)


def main():
    logging.getLogger("load_supplement").setLevel(logging.ERROR)
    cfg = load_config()
    xlsx = ROOT / cfg["paper"]["supplement_file"]
    freqs, fig1e, labels = ls.load_freq_sheet(xlsx, "Fig 1e")
    reported = ls.load_reported(xlsx)
    feats = pd.read_csv(ROOT / "data" / "features.csv").set_index("label").loc[labels]
    fits = pd.read_csv(ROOT / "data" / "fig1e_fits.csv").set_index("label").loc[labels]

    i1, i15 = int(np.flatnonzero(freqs == 1.0)[0]), int(np.flatnonzero(freqs == 1.5)[0])
    flip = []
    for lab, y in zip(labels, fig1e):
        r = knee_curve(fits.loc[lab, ["b_fit", "k_fit", "chi_fit"]].to_numpy(float), freqs) - np.log10(y)
        flip.append(bool(np.sign(r[i1]) != np.sign(r[i15])))

    df = pd.DataFrame({"group": [lab[:2] for lab in labels],
                       "pd_refit": feats["periodic_delta_peak_knee"].to_numpy(float),
                       "pd_reported": reported.loc[labels, "periodic_delta"].to_numpy(float),
                       "knee_flag": feats["knee_flag"].astype(bool).to_numpy(),
                       "age_lt_3": (feats["age"] < 3).to_numpy(),
                       "unclassified": (fits["outcome"] == "unclassified").to_numpy(),
                       "flip": np.array(flip)}, index=labels)
    if set(fits["outcome"]) - {"B", "unclassified"}:
        raise RuntimeError("Fig 1e outcomes other than B/unclassified; the 'unclassified' contrast assumes none")

    gap_subsets(df).to_csv(ROOT / "tables" / "periodic_gap_subsets.csv", index=False, float_format="%.5g")
    crosstabs(df).to_csv(ROOT / "tables" / "fig1e_crosstabs.csv", index=False, float_format="%.5g")


if __name__ == "__main__":
    main()
