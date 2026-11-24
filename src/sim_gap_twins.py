"""Why the synthetic baseline overpredicted the TD fixed - knee gap (post hoc, exploratory).

Design and readings were declared in NOTES (2026-11-16) before this was run. The gap: synthetic
(a) at nlv 0.02 gives a paired fixed - knee periodic delta of +0.321 (median); real TD gives +0.110.

Step 1, no new fits: TD paired d = fixed - knee periodic delta (peak model) by knee-frequency
bin of the knee fit, and in recordings near the synthetic baseline (knee 3.0-4.2 Hz, chi 2.0-2.4).

Step 2, twins: each TD recording rebuilt from its own knee-mode fit, refit in both modes.
    M0   median (b0, k0, chi0), no peaks, white noise SD 0.02      (synthetic (a))
    T00  own aperiodic, no peaks, white noise
    T10  own aperiodic + own knee-fit peak model, white noise
    T01  own aperiodic + own knee-model residual (log10 PSD - aperiodic - peaks)
    T11  own aperiodic + peaks + residual: the real spectrum (closure check)
White noise for recording j (features order) uses seed random_seed + j, the noise vectors of
(a)'s replicates 0-184 at nlv 0.02 (src/simulate.py). Effects on median d add up to M0 - T11:
heterogeneity M0 - T00; peaks, mean of (T00 - T10) and (T01 - T11); residual, mean of
(T00 - T01) and (T10 - T11). Interaction (T00 - T10) - (T01 - T11) is reported separately.

Follow-up (declared in NOTES 2026-11-16 before its run): the residual effect split into size and
structure with two more cells, white noise rescaled to the recording's own residual SD:
    T0s  own aperiodic, no peaks, T00's noise vector rescaled to SD = SD(own residual)
    T1s  own aperiodic + own peaks, the same noise
size = mean of (T00 - T0s) and (T10 - T1s); structure = mean of (T0s - T01) and (T1s - T11).
They add up to the residual effect. "Structure" is everything about the real residual except
its overall SD (correlation between bins, distribution shape, where in frequency it is large).
The twins inherit the knee-mode decomposition: they explain the generator's miss, not which
aperiodic model is right.

Reads data/features.csv and data/fits.npz (src/fit_spectra.py). Writes data/sim_gap_twins.csv
(per recording, local), tables/sim_gap_twins.csv and tables/sim_gap_by_knee_freq.csv (aggregates).
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from analysis import _pct, _settings
from config import load_config
from fit_spectra import fit_one
from simulate import generate

log = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
NLV = 0.02
CELLS = ("M0", "T00", "T10", "T01", "T11", "T0s", "T1s")
KNEE_BINS = [0, 2, 3, 4.2, 6, np.inf]
NEAR_BASELINE = {"knee_freq": (3.0, 4.2), "exponent": (2.0, 2.4)}
KNEE_RISE = 0.05  # split threshold for the descriptive T00 -> T11 rows, chosen after the cells were seen
LABEL = "post hoc, exploratory (NOTES 2026-11-16)"


def twin_cells(ap, peak, resid, ap_median, noise):
    """log10 spectra of the seven cells for one recording (all arrays over the same bins)."""
    noise_own = noise * (resid.std() / noise.std())  # same vector, sample SD = residual SD
    return {"M0": ap_median + noise, "T00": ap + noise, "T10": ap + peak + noise,
            "T01": ap + resid, "T11": ap + peak + resid,
            "T0s": ap + noise_own, "T1s": ap + peak + noise_own}


def effects(m):
    """Decomposition of M0 - T11 from per-cell values m (scalars or arrays over resamples)."""
    peaks = ((m["T00"] - m["T10"]) + (m["T01"] - m["T11"])) / 2
    resid = ((m["T00"] - m["T01"]) + (m["T10"] - m["T11"])) / 2
    size_no_peaks, size_with_peaks = m["T00"] - m["T0s"], m["T10"] - m["T1s"]
    return {"total": m["M0"] - m["T11"], "heterogeneity": m["M0"] - m["T00"], "peaks": peaks,
            "residual": resid, "interaction": (m["T00"] - m["T10"]) - (m["T01"] - m["T11"]),
            "residual_size": (size_no_peaks + size_with_peaks) / 2,
            "residual_structure": ((m["T0s"] - m["T01"]) + (m["T1s"] - m["T11"])) / 2,
            "residual_size_no_peaks": size_no_peaks, "residual_size_with_peaks": size_with_peaks}


def fit_twins(td, z, cfg):
    """Fit every cell of every TD recording in both modes; one row per (recording, cell)."""
    freqs, labels = z["freqs"], list(z["labels"])
    seed0 = int(cfg["analysis"]["random_seed"])
    b0, k0, chi0 = (float(td[c].median()) for c in ("offset_knee", "knee_knee", "exponent_knee"))
    ap_median = np.log10(generate(b0, k0, chi0, freqs=freqs))
    rows = []
    for j, lab in enumerate(td["label"]):
        i = labels.index(lab)
        ap, peak, log_psd = z["aperiodic_knee"][i], z["peak_knee"][i], z["log_psd"][i]
        r = td.iloc[j]
        own = np.log10(generate(r["offset_knee"], r["knee_knee"], r["exponent_knee"], freqs=freqs))
        if np.abs(own - ap).max() > 1e-6:
            raise RuntimeError(f"{lab}: stored knee aperiodic != knee curve from features.csv")
        noise = np.random.RandomState(seed0 + j).normal(0, NLV, len(freqs))
        for cell, log_spec in twin_cells(ap, peak, log_psd - ap - peak, ap_median, noise).items():
            row = {"label": lab, "cell": cell, "knee_freq": r["knee_freq_knee"], "exponent": r["exponent_knee"],
                   "resid_sd": float((log_psd - ap - peak).std())}
            for mode in ("knee", "fixed"):
                res = fit_one(freqs, 10 ** log_spec, mode, cfg)
                row.update({f"success_{mode}": res["success"],
                            f"pd_{mode}": res["periodic_delta_peak"],
                            f"n_peaks_below_6hz_{mode}": res["n_peaks_below_6hz"],
                            f"exponent_{mode}": res["exponent"], f"error_{mode}": res["error"]})
            rows.append(row)
        if (j + 1) % 25 == 0:
            log.info("%d/%d recordings", j + 1, len(td))
    out = pd.DataFrame(rows)
    out["d"] = out["pd_fixed"] - out["pd_knee"]
    return out, (b0, k0, chi0)


def cell_table(tw, cfg):
    """Median d (and fixed, knee) per cell, plus the effects, with paired bootstrap CIs."""
    n_boot, seed, level = _settings(cfg)
    wide = {c: tw[tw["cell"] == c].set_index("label") for c in CELLS}
    order = wide["T11"].index
    n = len(order)
    idx = np.random.default_rng(seed).integers(0, n, (n_boot, n))
    rows, boot = [], {}
    t00 = wide["T00"].loc[order]
    for c in CELLS:
        w = wide[c].loc[order]
        d = w["d"].to_numpy(float)
        boot[c] = np.median(d[idx], axis=1)
        lo, hi = _pct(boot[c], level)
        row = {"kind": "cell", "name": c, "n": n, "median_d": np.median(d), "ci_lower": lo, "ci_upper": hi,
               "frac_d_gt_0": float((d > 0).mean()),
               "median_pd_fixed": w["pd_fixed"].median(), "median_pd_knee": w["pd_knee"].median()}
        for m in ("fixed", "knee"):  # spread, because the median d and the difference of medians come apart
            for q in (0.25, 0.75, 0.9):
                row[f"pd_{m}_q{int(q * 100)}"] = w[f"pd_{m}"].quantile(q)
            delta = w[f"pd_{m}"] - t00[f"pd_{m}"]  # paired change against T00; d changes by dfixed - dknee
            row[f"median_d{m}_vs_T00"], row[f"mean_d{m}_vs_T00"] = delta.median(), delta.mean()
        row.update({"median_error_knee": w["error_knee"].median(), "median_error_fixed": w["error_fixed"].median(),
                    "median_n_peaks_below_6hz_fixed": w["n_peaks_below_6hz_fixed"].median(),
                    "median_n_peaks_below_6hz_knee": w["n_peaks_below_6hz_knee"].median(),
                    "success_both": float((w["success_knee"] & w["success_fixed"]).mean())})
        rows.append(row)
    point = effects({c: r["median_d"] for c, r in zip(CELLS, rows)})
    eboot = effects(boot)
    for name, v in point.items():
        lo, hi = _pct(eboot[name], level)
        rows.append({"kind": "effect", "name": name, "n": n, "median_d": v, "ci_lower": lo, "ci_upper": hi,
                     "share_of_total": v / point["total"],
                     "share_of_residual": v / point["residual"] if name.startswith("residual_") else np.nan})
    # T00 -> T11 split by whether knee mode gains periodic delta (threshold chosen after seeing the cells)
    real = wide["T11"].loc[order]
    dknee, dfixed = real["pd_knee"] - t00["pd_knee"], real["pd_fixed"] - t00["pd_fixed"]
    for name, m in ((f"T11_vs_T00_knee_rose_gt_{KNEE_RISE}", dknee > KNEE_RISE),
                    (f"T11_vs_T00_knee_rose_le_{KNEE_RISE}", dknee <= KNEE_RISE)):
        rows.append({"kind": "split", "name": name, "n": int(m.sum()), "median_d": real.loc[m, "d"].median(),
                     "median_d_T00": t00.loc[m, "d"].median(),
                     "median_dknee_vs_T00": dknee[m].median(), "median_dfixed_vs_T00": dfixed[m].median(),
                     "median_dd_vs_T00": (real.loc[m, "d"] - t00.loc[m, "d"]).median()})
    out = pd.DataFrame(rows)
    out["label"] = LABEL
    return out


def by_knee_freq(td, tw, cfg):
    """Step 1 (real d by knee-frequency bin and near the baseline), with the T00 twin alongside."""
    n_boot, seed, level = _settings(cfg)
    real = td.set_index("label")
    t00 = tw[tw["cell"] == "T00"].set_index("label").loc[real.index]
    kf = real["knee_freq_knee"]
    bins = pd.cut(kf, KNEE_BINS, right=False).astype(str).where(kf.notna(), "no knee frequency")
    near = (kf.between(*NEAR_BASELINE["knee_freq"]) & real["exponent_knee"].between(*NEAR_BASELINE["exponent"]))
    groups = [(b, bins == b) for b in [str(x) for x in pd.cut(kf, KNEE_BINS, right=False).cat.categories]]
    groups += [("no knee frequency", bins == "no knee frequency"), ("near baseline: knee 3.0-4.2 Hz, chi 2.0-2.4", near),
               ("all", pd.Series(True, index=real.index))]
    rows = []
    for name, m in groups:
        d = (real.loc[m, "periodic_delta_peak_fixed"] - real.loc[m, "periodic_delta_peak_knee"]).to_numpy(float)
        rng_idx = np.random.default_rng(seed).integers(0, len(d), (n_boot, len(d)))
        lo, hi = _pct(np.median(d[rng_idx], axis=1), level)
        rows.append({"subset": name, "n": int(m.sum()), "median_d": np.median(d), "ci_lower": lo, "ci_upper": hi,
                     "median_pd_fixed": real.loc[m, "periodic_delta_peak_fixed"].median(),
                     "median_pd_knee": real.loc[m, "periodic_delta_peak_knee"].median(),
                     "median_exponent_knee": real.loc[m, "exponent_knee"].median(),
                     "t00_median_d": t00.loc[m, "d"].median()})
    out = pd.DataFrame(rows)
    out["label"] = LABEL
    return out


def main():
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    cfg = load_config()
    feats = pd.read_csv(ROOT / "data" / "features.csv")
    td = feats[feats["group"] == "TD"].reset_index(drop=True)
    z = np.load(ROOT / "data" / "fits.npz")

    tw, (b0, k0, chi0) = fit_twins(td, z, cfg)
    log.info("M0 at b0=%.6f k0=%.6f chi0=%.6f (knee %.3f Hz)", b0, k0, chi0, k0 ** (1 / chi0))

    real_d = (td["periodic_delta_peak_fixed"] - td["periodic_delta_peak_knee"]).to_numpy(float)
    t11 = tw[tw["cell"] == "T11"].set_index("label").loc[td["label"], "d"].to_numpy(float)
    log.info("closure: T11 vs real d, max |diff| %.2e; medians %.6f vs %.6f",
             np.abs(t11 - real_d).max(), np.median(t11), np.median(real_d))
    if abs(np.median(t11) - np.median(real_d)) > 1e-4:
        raise RuntimeError("T11 (rebuilt real spectrum) does not reproduce the real TD median d")

    tw.to_csv(ROOT / "data" / "sim_gap_twins.csv", index=False)
    cell_table(tw, cfg).to_csv(ROOT / "tables" / "sim_gap_twins.csv", index=False, float_format="%.5g")
    by_knee_freq(td, tw, cfg).to_csv(ROOT / "tables" / "sim_gap_by_knee_freq.csv", index=False, float_format="%.5g")
    log.info("wrote data/sim_gap_twins.csv (local), tables/sim_gap_twins.csv, tables/sim_gap_by_knee_freq.csv")


if __name__ == "__main__":
    main()
