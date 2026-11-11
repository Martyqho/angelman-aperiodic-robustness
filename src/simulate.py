"""Synthetic validation: does the paper's periodic quantity separate peaks from aperiodic shifts?

Generator (matches specparam's knee model exactly):
    log10 P(f) = b - log10(k + f^chi) + G(f) + eps(f)
    G(f) = A * exp(-(f - cf)^2 / (2 sigma^2))     sigma is a standard deviation;
                                                   specparam reports bandwidth = 2 sigma
    eps ~ Normal(0, nlv) in log10 power            (white; not the chi-square noise of Welch)
on the data's grid, 1.0-49.5 Hz at 0.5 Hz (98 bins).

Baseline (b0, k0, chi0) = medians of knee-mode fit_one() over the 185 TD spectra.
Shifts come from the reported values: db and dchi are AS - TD medians of the Fig 2
offset and exponent; the peak height A makes the injected peak's band-rule value
equal the AS - TD median of the reported periodic delta. A/2 and 2A give a range.

Conditions: (a) baseline; (b) offset only; (c) exponent only (k fixed, so the knee
frequency k^(1/chi) moves down); (c') exponent only with k = k0^(chi1/chi0), so the
knee frequency stays at k0^(1/chi0); (d) peak only; (e) all three. Every spectrum is
fit in knee and fixed mode with fit_one().

Seeds are paired across conditions: replicate r at noise level i uses
seed = random_seed + 1000 * i + r, so (b) - (a) etc. differ only by the condition.

Writes data/sim_spectra.csv (per spectrum, local), tables/sim_summary.csv and
figures/fig_synthetic_validation.png.

Usage: python src/simulate.py [n_replicates]
"""

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import load_supplement as ls
from config import load_config
from fit_spectra import DELTA, band_power, fit_one

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
XLSX = ROOT / "refs" / "supplementary-data-1.xlsx"
DATA, TABLES, FIGURES = ROOT / "data", ROOT / "tables", ROOT / "figures"

FREQS = np.arange(1.0, 49.5 + 0.25, 0.5)
PEAK_CF, PEAK_SIGMA = 3.0, 0.5
NOISE_LEVELS = (0.02, 0.05)
A_SCALES = (0.5, 1.0, 2.0)
N_REPLICATES = 200
MODES = ("knee", "fixed")
PEAK_HIT_HZ = 0.5  # a fitted peak within +-0.5 Hz of the injected 3 Hz counts as a hit
METRICS = ["raw_delta", "periodic_delta_peak", "periodic_delta_resid", "offset", "exponent",
           "knee_freq", "r2", "error", "n_peaks", "n_peaks_below_6hz",
           "offset_err", "exponent_err", "knee_freq_err"]


def generate(b, k, chi, A=0.0, nlv=0.0, seed=0, freqs=FREQS, cf=PEAK_CF, sigma=PEAK_SIGMA):
    """Linear PSD from the knee model plus one Gaussian peak plus white log10 noise."""
    log_p = b - np.log10(k + freqs ** chi) + A * np.exp(-(freqs - cf) ** 2 / (2 * sigma ** 2))
    if nlv > 0:
        log_p = log_p + np.random.RandomState(seed).normal(0, nlv, len(freqs))
    return 10 ** log_p


def peak_band_value(A, freqs=FREQS):
    """Band-rule 2-4 Hz value of the injected peak alone (log10 units)."""
    return float(band_power(freqs, A * np.exp(-(freqs - PEAK_CF) ** 2 / (2 * PEAK_SIGMA ** 2)), *DELTA))


def baseline_from_td(cfg, path=XLSX):
    """Medians of knee-mode fits over the 185 TD spectra."""
    freqs, psd, labels = ls.load_spectra(path)
    td = [i for i, lab in enumerate(labels) if lab.startswith("TD")]
    fits = pd.DataFrame([fit_one(freqs, psd[i], "knee", cfg) for i in td])
    med = fits[["offset", "knee", "exponent", "knee_freq"]].median()
    log.info("TD knee fits: %d/%d succeeded, %d with knee <= 0 (no knee frequency)",
             fits["success"].sum(), len(td), (fits["knee"] <= 0).sum())
    log.info("baseline b0=%.4f k0=%.4f chi0=%.4f; median knee freq %.3f Hz (k0^(1/chi0) = %.3f Hz)",
             med["offset"], med["knee"], med["exponent"], med["knee_freq"],
             med["knee"] ** (1 / med["exponent"]))
    return {"b0": med["offset"], "k0": med["knee"], "chi0": med["exponent"],
            "knee_freq_median": med["knee_freq"], "n_td": len(td),
            "n_positive_knee": int((fits["knee"] > 0).sum())}


def shifts_from_reported(path=XLSX):
    """db, dchi from Fig 2 (AS - NT medians); A from the reported periodic delta."""
    table = ls.load_table(path)
    med = table.groupby("genotype")[["offset", "exponent"]].median()
    reported = ls.load_reported(path)
    pd_med = reported.groupby(reported.index.str[:2])["periodic_delta"].median()
    d_pd = pd_med["AS"] - pd_med["TD"]
    A = d_pd / peak_band_value(1.0)
    out = {"db": med.loc["AS", "offset"] - med.loc["NT", "offset"],
           "dchi": med.loc["AS", "exponent"] - med.loc["NT", "exponent"],
           "d_periodic_delta": d_pd, "A": A}
    log.info("shifts db=%.4f dchi=%.4f; periodic delta AS-TD=%.4f -> A=%.4f",
             out["db"], out["dchi"], d_pd, A)
    return out


def conditions(base, shifts):
    """[(condition, a_scale, params)] for the grid."""
    b0, k0, chi0 = base["b0"], base["k0"], base["chi0"]
    b1, chi1, A = b0 + shifts["db"], chi0 + shifts["dchi"], shifts["A"]
    out = [("a_baseline", 0.0, dict(b=b0, k=k0, chi=chi0, A=0.0)),
           ("b_offset", 0.0, dict(b=b1, k=k0, chi=chi0, A=0.0)),
           ("c_exponent", 0.0, dict(b=b0, k=k0, chi=chi1, A=0.0)),
           ("c_exponent_kf", 0.0, dict(b=b0, k=k0 ** (chi1 / chi0), chi=chi1, A=0.0))]
    for s in A_SCALES:
        out.append(("d_peak", s, dict(b=b0, k=k0, chi=chi0, A=s * A)))
    for s in A_SCALES:
        out.append(("e_all", s, dict(b=b1, k=k0, chi=chi1, A=s * A)))
    return out


def run_grid(cfg, base, shifts, n_rep=N_REPLICATES):
    seed0 = int(cfg["analysis"]["random_seed"])
    rows = []
    for i, nlv in enumerate(NOISE_LEVELS):
        for cond, a_scale, p in conditions(base, shifts):
            truth_raw = float(band_power(FREQS, generate(**p), *DELTA))
            truth_pd = peak_band_value(p["A"])
            for r in range(n_rep):
                seed = seed0 + 1000 * i + r
                psd = generate(**p, nlv=nlv, seed=seed)
                for mode in MODES:
                    res = fit_one(FREQS, psd, mode, cfg)
                    cfs = np.array(res.pop("peak_cfs"))
                    near = cfs[np.argmin(np.abs(cfs - PEAK_CF))] if len(cfs) else np.nan
                    rows.append({"condition": cond, "a_scale": a_scale, "nlv": nlv, "replicate": r,
                                 "seed": seed, "mode": mode, "true_b": p["b"], "true_k": p["k"], "true_chi": p["chi"],
                                 "true_knee_freq": p["k"] ** (1 / p["chi"]),
                                 "true_A": p["A"], "true_raw_delta": truth_raw,
                                 "true_periodic_delta": truth_pd, **res,
                                 "peak_cf_nearest_3hz": near,
                                 "peak_hit_3hz": bool(np.abs(near - PEAK_CF) <= PEAK_HIT_HZ)})
        log.info("nlv=%.2f done (%d fits so far)", nlv, len(rows))
    return pd.DataFrame(rows)


def _q(x, q):
    return float(np.nanpercentile(x, q)) if np.isfinite(x).any() else float("nan")


def _truth(g):
    """Generating parameters of one grid cell, so the summary table carries its provenance."""
    first = g.iloc[0]
    return {"true_b": first["true_b"], "true_k": first["true_k"], "true_chi": first["true_chi"],
            "true_knee_freq": first["true_k"] ** (1 / first["true_chi"]), "true_A": first["true_A"]}


def summarize(sims):
    """One row per (condition, a_scale, nlv, mode), plus 'fixed-knee' paired differences.

    Parameter recovery: *_err = estimated - true (offset vs b, exponent vs chi, knee
    frequency vs k^(1/chi); knee frequency is NaN in fixed mode, and the fixed-mode
    offset is a different parameter from b). Peak detection: frac_any_peak_below_6hz
    (a false-positive rate in the no-peak cells) and frac_peak_hit_3hz (a fitted peak
    within +-0.5 Hz of 3 Hz; recovery in (d)/(e), a false positive elsewhere).
    """
    keys = ["condition", "a_scale", "nlv"]
    sims = sims.assign(offset_err=sims["offset"] - sims["true_b"],
                       exponent_err=sims["exponent"] - sims["true_chi"],
                       knee_freq_err=sims["knee_freq"] - sims["true_k"] ** (1 / sims["true_chi"]))
    base = sims[sims["condition"] == "a_baseline"].groupby(["nlv", "mode"])["periodic_delta_peak"]
    lo, hi = base.quantile(0.025), base.quantile(0.975)

    rows = []
    for (cond, a_scale, nlv, mode), g in sims.groupby(keys + ["mode"]):
        row = {"condition": cond, "a_scale": a_scale, "nlv": nlv, "mode": mode, **_truth(g), "n": len(g),
               "success_rate": g["success"].mean(),
               "true_raw_delta": g["true_raw_delta"].iloc[0],
               "true_periodic_delta": g["true_periodic_delta"].iloc[0]}
        for m in METRICS:
            row[f"{m}_median"] = _q(g[m], 50)
            row[f"{m}_q25"], row[f"{m}_q75"] = _q(g[m], 25), _q(g[m], 75)
        row["frac_any_peak_below_6hz"] = float((g["n_peaks_below_6hz"] >= 1).mean())
        if "peak_hit_3hz" in g:
            row["frac_peak_hit_3hz"] = float(g["peak_hit_3hz"].astype(bool).mean())
        pdp = g["periodic_delta_peak"]
        row["frac_pd_peak_in_baseline_95"] = float(((pdp >= lo[(nlv, mode)]) & (pdp <= hi[(nlv, mode)])).mean())
        rows.append(row)

    wide = sims.pivot_table(index=keys + ["replicate"], columns="mode",
                            values=["periodic_delta_peak", "periodic_delta_resid", "exponent"])
    truth = sims.groupby(keys)
    for (cond, a_scale, nlv), g in wide.groupby(level=keys):
        row = {"condition": cond, "a_scale": a_scale, "nlv": nlv, "mode": "fixed-knee",
               **_truth(truth.get_group((cond, a_scale, nlv))), "n": len(g)}
        for m in ["periodic_delta_peak", "periodic_delta_resid", "exponent"]:
            d = (g[(m, "fixed")] - g[(m, "knee")]).to_numpy()
            row[f"{m}_median"], row[f"{m}_q25"], row[f"{m}_q75"] = _q(d, 50), _q(d, 25), _q(d, 75)
        rows.append(row)
    return pd.DataFrame(rows)


def main(n_rep=N_REPLICATES):
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    logging.getLogger("load_supplement").setLevel(logging.ERROR)
    cfg = load_config()
    base = baseline_from_td(cfg)
    shifts = shifts_from_reported()
    sims = run_grid(cfg, base, shifts, n_rep)
    summary = summarize(sims)

    for d in (DATA, TABLES, FIGURES):
        d.mkdir(exist_ok=True)
    sims.to_csv(DATA / "sim_spectra.csv", index=False)
    summary.to_csv(TABLES / "sim_summary.csv", index=False, float_format="%.5g")
    log.info("wrote %s (local), %s", DATA / "sim_spectra.csv", TABLES / "sim_summary.csv")

    from plot_simulation import plot
    plot(sims, base, shifts, FIGURES / "fig_synthetic_validation.png")
    log.info("wrote %s", FIGURES / "fig_synthetic_validation.png")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else N_REPLICATES)
