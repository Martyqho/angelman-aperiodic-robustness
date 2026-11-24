"""Secondary (a): agreement between the knee refit and the paper's reported values.

Reads data/features.csv (python src/fit_spectra.py) and the reported values from
Supplementary Data 1 ('Fig 1b,d,f,g'). For exponent, offset and periodic delta
(refit: band rule over the knee-mode fitted peak model), per group (TD, AS):
n, Pearson r, mean bias, SD, 95% limits of agreement (bias +- 1.96 SD), median
of refit - reported, median |refit - reported| and the fraction with |difference| > 0.15. 'all' includes every converged recording (the analysis set);
'no_knee_flag' excludes knee-flagged recordings (sensitivity).

Writes tables/validation.csv (aggregate) and tables/deviations.csv, which lists two things:
- every recorded difference between this project's values and the paper's (kind 'reported
  value': AS159 Total 1-50 Hz, and the knee-refit biases for exponent, offset and periodic
  delta, analysis set 'all'), computed here;
- every difference from PROTOCOL v2.3 and why (all other kinds: procedure, specification,
  amendment, expectation, addition, exploratory plan). These rows are kept by hand in
  PROTOCOL_DEVIATIONS below; each 'why' is what NOTES recorded at the time, or says that no
  reason was recorded.
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


# Differences from PROTOCOL v2.3, written at the Phase 4 close (2026-11-16). Text only.
PROTOCOL_DEVIATIONS = [
    {"kind": "procedure", "item": "Primary: paired plots", "group": "all",
     "protocol_says": "median paired difference and bootstrap CI, plus paired plots",
     "what_was_done": "The committed figure bins recordings by pair mean (quintiles per group) and shows the "
                      "median and IQR of d with the group median line and CI band. The per-recording paired "
                      "plot is local (data/fig_primary_paired_points.png).",
     "why": "A per-recording plot in a committed file would redistribute derived per-recording values "
            "(CC BY-NC-ND 4.0).",
     "source": "NOTES 2026-11-14 (Figure); figures/fig_primary_paired.png"},
    {"kind": "procedure", "item": "Secondary (b): where it is reported", "group": "TD, AS",
     "protocol_says": "(b) AS-vs-TD Cohen's d for the primary periodic quantity and for exponent under each mode",
     "what_was_done": "Reported in tables/group_contrast_by_mode.csv, labelled predeclared. tables/secondary_b.csv "
                      "holds the fixed - knee change in exponent and offset by group, labelled exploratory.",
     "why": "The Phase 4 analysis spec defined (b) as the fixed - knee parameter changes; both were kept and "
            "every row is labelled.",
     "source": "NOTES 2026-11-14 (Spec vs PROTOCOL)"},
    {"kind": "procedure", "item": "Secondary (c): GAM shape comparison", "group": "TD",
     "protocol_says": "linear first; a GAM only to compare shape with the paper's Fig. 3 TD curve",
     "what_was_done": "Not run. Only linear slopes on log age were fitted, so whether the paper's TD curve keeps "
                      "its shape under fixed mode is answered for the linear slope only.",
     "why": "No reason recorded (NOTES 2026-11-14 says only 'No GAM').",
     "source": "NOTES 2026-11-14; tables/secondary_c.csv"},
    {"kind": "specification", "item": "Statistics: log age", "group": "TD (AS exploratory)",
     "protocol_says": "linear model of each feature on log age per mode",
     "what_was_done": "log2(age); slopes are per doubling of age.",
     "why": "PROTOCOL does not fix the base of the logarithm; recorded with the choice.",
     "source": "NOTES 2026-11-14; tables/secondary_c.csv"},
    {"kind": "amendment", "item": "QC gate item 2", "group": "all",
     "protocol_says": "v2.2 wording: both fits succeed for at least 95% of recordings",
     "what_was_done": "Failed as written on 2026-11-04 when knee flags were counted as failures (knee mode "
                      "319/344 = 92.7%). Resolved 2026-11-05: item 2 means convergence (344/344 in both modes); "
                      "wording corrected in v2.3.",
     "why": "Excluding knee-flagged recordings conditions on the knee-mode outcome, one arm of the comparison "
            "(decision recorded 2026-11-05).",
     "source": "NOTES 2026-11-04 (gate verdict), 2026-11-05 (Item 2 resolved); PROTOCOL change log v2.2-v2.3"},
    {"kind": "expectation", "item": "Expectation 2026-11-03: size in TD", "group": "TD",
     "protocol_says": "fixed-mode periodic delta exceeds knee-mode by roughly 0.3 log10 units in TD",
     "what_was_done": "+0.110 [0.094, 0.129]. Direction held; size not held (the CI excludes 0.3).",
     "why": "Post hoc (NOTES 2026-11-16): the synthetic generator lacked each recording's own peaks (45% of the "
            "0.21 gap) and the structure of the real residual (32%); residual size 8%, the median-curve "
            "shortcut 15%.",
     "source": "tables/primary.csv; tables/sim_gap_twins.csv; NOTES 2026-11-14, 2026-11-16"},
    {"kind": "expectation", "item": "Expectation 2026-11-03: mechanism", "group": "AS vs TD",
     "protocol_says": "because AS knee frequencies are lower and exponents steeper; if the data do that, "
                      "the result is a confirmed mechanism",
     "what_was_done": "Not tested. The predicted direction held (larger gap in AS; the AS-TD contrast increased "
                      "under fixed mode), so the mechanism is neither confirmed nor ruled out.",
     "why": "The primary and group-contrast numbers measure the change in the contrast, not its cause.",
     "source": "NOTES 2026-11-14 (Expectations, clause by clause)"},
    {"kind": "expectation", "item": "Expectation 2026-11-05: flagged recordings and the age trend", "group": "TD",
     "protocol_says": "knee-flagged recordings may themselves produce an age trend in the fixed - knee exponent "
                      "difference",
     "what_was_done": "Not supported. TD slope -0.241 per doubling of age with them, -0.227 without (94%).",
     "why": "The trend is in the unflagged recordings.",
     "source": "tables/secondary_c.csv; NOTES 2026-11-14"},
    {"kind": "addition", "item": "Exploratory additions in Phase 4", "group": "TD, AS",
     "protocol_says": "not in PROTOCOL",
     "what_was_done": "AS rows of (c); AS - TD median periodic delta per mode and its fixed - knee change; "
                      "fixed - knee change in exponent and offset, with knee-flagged and unflagged subsets; "
                      "mean d with CI.",
     "why": "Requested in the Phase 4 analysis spec (parameter changes; mean d, labelled descriptive; the AS - TD "
            "median difference, for the email and README) or added to test the 2026-11-05 expectation (flag "
            "subsets). AS slopes: computed alongside TD, no separate reason recorded. Every row is labelled "
            "exploratory or descriptive.",
     "source": "tables/secondary_b.csv, secondary_c.csv, group_contrast_by_mode.csv, primary.csv; NOTES 2026-11-14"},
    {"kind": "addition", "item": "Post hoc follow-ups", "group": "TD, AS",
     "protocol_says": "not in PROTOCOL",
     "what_was_done": "Gap decomposition in (d), periodic-gap subsets and Fig 1e cross-tabs (2026-11-14); "
                      "synthetic-vs-real TD gap twins with the residual size/structure split (2026-11-16).",
     "why": "Requested after the results they follow up; labelled post hoc. Readings for the twins were "
            "declared in NOTES before they ran.",
     "source": "tables/secondary_d.csv, periodic_gap_subsets.csv, fig1e_crosstabs.csv, sim_gap_twins.csv, "
               "sim_gap_by_knee_freq.csv; NOTES 2026-11-14, 2026-11-16"},
    {"kind": "exploratory plan", "item": "PROTOCOL's exploratory list", "group": "all",
     "protocol_says": "Exploratory, not predeclared: fixed mode 2-40 Hz; doublexp; the Fig 6 mouse PCA; HBN",
     "what_was_done": "None run.",
     "why": "No reason recorded.",
     "source": "PROTOCOL v2.3 (Exploratory, not predeclared)"},
]
TEXT_COLUMNS = ["kind", "protocol_says", "what_was_done", "why"]


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
    rows = [{"kind": "reported value", "item": "Total 1-50 Hz", "group": "AS159", "measure": "value (0.5 x sum of 98 bins)",
             "reported": reported.loc["AS159", "total"], "ours": ours_total,
             "difference": ours_total - reported.loc["AS159", "total"],
             "relative": ours_total / reported.loc["AS159", "total"] - 1,
             "source": "src/validate.py check 1; only recording of 344 off by > 1e-5"}]
    v = validation[validation["subset"] == "all"]
    for _, r in v.iterrows():
        rows.append({"kind": "reported value", "item": f"{r['feature']}, knee refit", "group": r["group"],
                     "measure": "mean over recordings; difference = refit - reported",
                     "reported": r["mean_reported"], "ours": r["mean_refit"],
                     "difference": r["mean_refit_minus_reported"], "relative": float("nan"),
                     "source": f"tables/validation.csv (n {r['n']}, LoA {r['loa_lower']:.3f} to "
                               f"{r['loa_upper']:.3f}, median |d| {r['median_abs_diff']:.3f})"})
    out = pd.DataFrame(rows + PROTOCOL_DEVIATIONS)
    first = ["item", "group", "measure", "reported", "ours", "difference", "relative", "source"]
    return out[first + TEXT_COLUMNS]


if __name__ == "__main__":
    main()
