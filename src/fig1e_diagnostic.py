"""Is Fig 1e one specparam knee fit? PROTOCOL secondary (a), continued (v2.3).

Per recording, fit b - log10(k + f^chi) to log10(Fig 1e) over the 98 bins with
scipy.optimize.least_squares, starting from the reported offset and exponent and
k = 3^chi (a 3 Hz knee). Outcomes (declared in NOTES.md 2026-11-05 before this ran):

  (A) max |residual| < 1e-4 and |b - offset|, |chi - exponent| both <= 1e-3
  (C) max |residual| < 1e-4 and either differs by more than 1e-3
  (B) max |residual| > 1e-3 and structured (lag-1 autocorrelation of the residual >= 0.5)
  unclassified: anything else

Side checks declared with it:
  1. reported periodic delta of TD150 and AS151 (printed)
  2. per group, does knee refit - reported periodic delta scale with reported periodic
     delta? OLS slope with 95% CI and Spearman rho, against reported and against
     the pair mean (refit + reported) / 2.

Writes data/fig1e_fits.csv (per recording, local), tables/fig1e_diagnostic.csv and
tables/periodic_gap_scaling.csv (aggregate). Needs data/features.csv (src/fit_spectra.py).
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import least_squares

import load_supplement as ls
from config import load_config

ROOT = Path(__file__).resolve().parent.parent
RESID_SINGLE, RESID_MULTI, PARAM_TOL, AC_STRUCTURED = 1e-4, 1e-3, 1e-3, 0.5
KNEE_GUESS_HZ = 3.0


def knee_curve(params, f):
    b, k, chi = params
    return b - np.log10(k + f ** chi)


def fit_curve(f, y, b0, chi0):
    x0 = [b0, KNEE_GUESS_HZ ** chi0, chi0]

    def resid(p):
        inner = p[1] + f ** p[2]
        return np.where(inner > 0, p[0] - np.log10(np.where(inner > 0, inner, 1.0)) - y, 1e3)

    sol = least_squares(resid, x0, method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    r = knee_curve(sol.x, f) - y
    return sol, r


def lag1_autocorr(r):
    r = r - r.mean()
    den = (r ** 2).sum()
    return float((r[:-1] * r[1:]).sum() / den) if den > 0 else float("nan")


def classify(max_r, ac, db, dchi):
    if max_r < RESID_SINGLE:
        return "A" if (abs(db) <= PARAM_TOL and abs(dchi) <= PARAM_TOL) else "C"
    if max_r > RESID_MULTI and ac >= AC_STRUCTURED:
        return "B"
    return "unclassified"


def scaling(d, x):
    lr = stats.linregress(x, d)
    t = stats.t.ppf(0.975, len(x) - 2)
    rho, p = stats.spearmanr(x, d)
    return {"slope": lr.slope, "slope_ci_lower": lr.slope - t * lr.stderr,
            "slope_ci_upper": lr.slope + t * lr.stderr, "intercept": lr.intercept,
            "spearman_rho": rho, "spearman_p": p, "n": len(x)}


def main():
    logging.getLogger("load_supplement").setLevel(logging.ERROR)
    cfg = load_config()
    xlsx = ROOT / cfg["paper"]["supplement_file"]
    freqs, fig1e, labels = ls.load_freq_sheet(xlsx, "Fig 1e")
    reported = ls.load_reported(xlsx)
    if list(reported.index) != labels:
        raise RuntimeError("Fig 1e and reported labels differ")
    if not (fig1e > 0).all():
        raise ValueError("Fig 1e: non-positive values")

    rows = []
    for lab, y_lin in zip(labels, fig1e):
        y = np.log10(y_lin)
        b0, chi0 = reported.loc[lab, "offset"], reported.loc[lab, "exponent"]
        sol, r = fit_curve(freqs, y, b0, chi0)
        b, k, chi = sol.x
        max_r, ac = float(np.abs(r).max()), lag1_autocorr(r)
        db, dchi = b - b0, chi - chi0
        rows.append({"label": lab, "group": lab[:2], "success": bool(sol.success), "b_fit": b, "k_fit": k,
                     "chi_fit": chi, "knee_freq_fit": k ** (1 / chi) if k > 0 else np.nan,
                     "offset_reported": b0, "exponent_reported": chi0, "db": db, "dchi": dchi,
                     "max_abs_resid": max_r, "argmax_resid_hz": float(freqs[np.abs(r).argmax()]),
                     "rms_resid": float(np.sqrt((r ** 2).mean())), "lag1_autocorr": ac,
                     "outcome": classify(max_r, ac, db, dchi)})
    fits = pd.DataFrame(rows)

    summary = []
    for g in ("TD", "AS", "all"):
        sub = fits if g == "all" else fits[fits["group"] == g]
        row = {"group": g, "n": len(sub), "converged": int(sub["success"].sum())}
        for o in ("A", "B", "C", "unclassified"):
            row[f"n_{o}"] = int((sub["outcome"] == o).sum())
        for q in (0, 0.05, 0.25, 0.5, 0.75, 0.95, 1):
            row[f"max_resid_q{int(q * 100):02d}"] = sub["max_abs_resid"].quantile(q)
        row["median_lag1_autocorr"] = sub["lag1_autocorr"].median()
        row["median_abs_db"] = sub["db"].abs().median()
        row["median_abs_dchi"] = sub["dchi"].abs().median()
        summary.append(row)
    summary = pd.DataFrame(summary)

    feats = pd.read_csv(ROOT / "data" / "features.csv").set_index("label")
    refit = feats.loc[labels, "periodic_delta_peak_knee"].to_numpy(float)
    rep = reported["periodic_delta"].to_numpy(float)
    grp = np.array([lab[:2] for lab in labels])
    scale = []
    for g in ("TD", "AS"):
        m = grp == g
        d = refit[m] - rep[m]
        for axis, x in (("reported", rep[m]), ("pair_mean", (refit[m] + rep[m]) / 2)):
            scale.append({"group": g, "x_axis": axis, **scaling(d, x)})
    scale = pd.DataFrame(scale)

    (ROOT / "data").mkdir(exist_ok=True)
    fits.to_csv(ROOT / "data" / "fig1e_fits.csv", index=False, float_format="%.10g")
    summary.to_csv(ROOT / "tables" / "fig1e_diagnostic.csv", index=False, float_format="%.5g")
    scale.to_csv(ROOT / "tables" / "periodic_gap_scaling.csv", index=False, float_format="%.5g")

    print("Fig 1e outcomes: " + "; ".join(
        f"{r['group']} A {r['n_A']}, B {r['n_B']}, C {r['n_C']}, unclassified {r['n_unclassified']}"
        for _, r in summary[summary["group"] != "all"].iterrows()))
    for lab in ("TD150", "AS151"):
        print(f"reported periodic delta {lab}: {reported.loc[lab, 'periodic_delta']:.6g}")


if __name__ == "__main__":
    main()
