"""One specparam fit with the paper's settings, the paper's band rule, and the full fit pass.

fit_one() is the only place a spectrum is fit; the fit pass over the 344
recordings (main(), below) and the synthetic validation (src/simulate.py) both call it.

Fit pass (python src/fit_spectra.py): all 344 Fig 1a spectra in knee and fixed mode.
Writes, all local (per-recording derived values, CC BY-NC-ND):
  data/features.csv  one row per recording: label, group, id, age, sex, del_status,
                     every fit_one scalar with a _knee / _fixed suffix, knee_flag
  data/peaks.csv     one row per fitted peak: label, mode, cf, pw (log10), bw (Hz)
  data/fits.npz      freqs, labels, and for each mode the aperiodic fit, peak model and
                     raw residual (log10, shape (344, 98)), plus log10 PSD
Prints convergence and knee-flag counts by group, nothing else.

knee_flag (knee mode): knee <= 0 (knee frequency undefined) or knee frequency
outside [0.5, 20] Hz. specparam 2.0.0rc2 fits the knee unbounded (ap_bounds are
+-inf), so "knee at bound" cannot occur and is not a separate flag.
"""

import sys
from pathlib import Path

import numpy as np
from specparam import SpectralModel

DELTA = (2.0, 4.0)
LOW_PEAK_HZ = 6.0
KNEE_FREQ_RANGE = (0.5, 20.0)


def band_power(freqs, x, lo, hi):
    """(hi - lo) * mean of x over bins lo <= f <= hi: the paper's band rule.

    Works along the last axis, so x may be one spectrum or (n, n_freqs).
    """
    m = (freqs >= lo) & (freqs <= hi)
    return (hi - lo) * np.asarray(x)[..., m].mean(axis=-1)


def model_settings(cfg):
    sp = cfg["specparam"]
    return dict(
        peak_width_limits=tuple(sp["peak_width_limits_hz"]),
        max_n_peaks=sp["max_n_peaks"],
        min_peak_height=sp["min_peak_height"],
        peak_threshold=sp["peak_threshold"],
    )


def fit_one(freqs, psd_linear, mode, cfg, return_components=False):
    """Fit one linear PSD in `mode` ('knee' or 'fixed') with the paper's settings.

    Returns a dict: success, offset, knee (NaN in fixed mode), exponent, knee_freq
    (k^(1/exponent), Hz), r2, error, n_peaks, n_peaks_below_6hz, peak_cfs (tuple of
    fitted peak center frequencies, Hz), and three 2-4 Hz
    band-rule values: periodic_delta_peak (fitted peak model, log10; the paper's
    quantity), periodic_delta_resid (log10(psd) - log10(aperiodic fit)) and
    raw_delta (linear PSD).

    With return_components=True, returns (dict, components): components has
    'aperiodic', 'peak' and 'residual' (log10, over the fit range; NaN if the fit
    failed) and 'peak_params' (n_peaks, 3: cf, pw, bw).
    """
    freqs = np.asarray(freqs, float)
    psd_linear = np.asarray(psd_linear, float)
    fm = SpectralModel(aperiodic_mode=mode, verbose=False, **model_settings(cfg))
    fm.fit(freqs, psd_linear, cfg["specparam"]["fit_range_hz"])

    out = {"raw_delta": float(band_power(freqs, psd_linear, *DELTA))}
    if not fm.has_model:
        nan = float("nan")
        out.update(success=False, offset=nan, knee=nan, exponent=nan, knee_freq=nan, r2=nan,
                   error=nan, n_peaks=0, n_peaks_below_6hz=0, peak_cfs=(),
                   periodic_delta_peak=nan, periodic_delta_resid=nan)
        if return_components:
            blank = np.full(len(freqs), nan)
            return out, {"aperiodic": blank, "peak": blank, "residual": blank,
                         "peak_params": np.empty((0, 3))}
        return out

    f = fm.freqs  # the fit range; equals `freqs` when the fit range covers all bins
    ap_fit = fm.get_model("aperiodic", "log")
    peak_fit = fm.get_model("peak", "log")
    log_psd = fm.power_spectrum
    # Component checks: the residual is specparam's flattened spectrum, and the
    # two components add up to the full model (both in log10).
    if not np.allclose(log_psd - ap_fit, fm._spectrum_flat):
        raise RuntimeError("log10(psd) - aperiodic fit != flattened spectrum")
    if not np.allclose(ap_fit + peak_fit, fm.get_model("full", "log")):
        raise RuntimeError("aperiodic fit + peak fit != modeled spectrum")

    ap = fm.aperiodic_params_
    offset, exponent = ap[0], ap[-1]
    knee = ap[1] if mode == "knee" else float("nan")
    cfs = fm.peak_params_[:, 0] if len(fm.peak_params_) else np.array([])
    out.update(
        success=True,
        offset=float(offset),
        knee=float(knee),
        exponent=float(exponent),
        knee_freq=float(knee ** (1 / exponent)) if mode == "knee" and knee > 0 else float("nan"),
        r2=float(fm.r_squared_),
        error=float(fm.error_),
        n_peaks=int(len(cfs)),
        n_peaks_below_6hz=int((cfs < LOW_PEAK_HZ).sum()),
        peak_cfs=tuple(float(c) for c in cfs),
        periodic_delta_peak=float(band_power(f, peak_fit, *DELTA)),
        periodic_delta_resid=float(band_power(f, log_psd - ap_fit, *DELTA)),
    )
    if return_components:
        return out, {"aperiodic": ap_fit, "peak": peak_fit, "residual": log_psd - ap_fit,
                     "peak_params": np.asarray(fm.peak_params_).reshape(-1, 3)}
    return out


def knee_flag(knee, knee_freq):
    """True where the knee fit has no usable knee frequency (see module docstring)."""
    lo, hi = KNEE_FREQ_RANGE
    knee, knee_freq = np.asarray(knee, float), np.asarray(knee_freq, float)
    return ~(knee > 0) | ~((knee_freq >= lo) & (knee_freq <= hi))


SCALARS = ["success", "offset", "knee", "exponent", "knee_freq", "r2", "error", "n_peaks",
           "n_peaks_below_6hz", "periodic_delta_peak", "periodic_delta_resid"]


def main():
    import logging

    import pandas as pd

    import load_supplement as ls
    from config import load_config

    logging.getLogger("load_supplement").setLevel(logging.ERROR)
    root = Path(__file__).resolve().parent.parent
    cfg = load_config()
    xlsx = root / cfg["paper"]["supplement_file"]
    freqs, psd, labels = ls.load_spectra(xlsx)
    tidy, _ = ls.load_all(xlsx)
    if list(tidy["label"]) != labels:
        raise RuntimeError("load_all and load_spectra disagree on order")
    if (cfg["specparam"]["fit_range_hz"][0] > freqs[0]) or (cfg["specparam"]["fit_range_hz"][1] < freqs[-1]):
        raise RuntimeError("fit range does not cover all bins; components would not align with freqs")

    feats = tidy[["label", "group", "id", "age", "sex", "del_status"]].copy()
    arrays = {"freqs": freqs, "labels": np.array(labels), "log_psd": np.log10(psd)}
    peaks = []
    for mode in ("knee", "fixed"):
        rows, comps = [], {k: [] for k in ("aperiodic", "peak", "residual")}
        for lab, x in zip(labels, psd):
            res, c = fit_one(freqs, x, mode, cfg, return_components=True)
            res.pop("peak_cfs")
            rows.append(res)
            for k in comps:
                comps[k].append(c[k])
            for cf, pw, bw in c["peak_params"]:
                peaks.append({"label": lab, "mode": mode, "cf": cf, "pw": pw, "bw": bw})
        df = pd.DataFrame(rows)
        if mode == "knee":
            feats["raw_delta"] = df["raw_delta"].to_numpy()
        feats = feats.join(df[SCALARS].add_suffix(f"_{mode}"))
        for k, v in comps.items():
            arrays[f"{k}_{mode}"] = np.vstack(v)
    feats["knee_flag"] = knee_flag(feats["knee_knee"], feats["knee_freq_knee"])

    data = root / "data"
    data.mkdir(exist_ok=True)
    feats.to_csv(data / "features.csv", index=False, float_format="%.10g")
    pd.DataFrame(peaks).to_csv(data / "peaks.csv", index=False, float_format="%.10g")
    np.savez_compressed(data / "fits.npz", **arrays)

    for mode in ("knee", "fixed"):
        ok = feats.groupby("group")[f"success_{mode}"].agg(["sum", "size"])
        print(f"converged, {mode}: " + ", ".join(f"{g} {int(r['sum'])}/{int(r['size'])}"
                                                for g, r in ok.loc[["TD", "AS"]].iterrows()))
    flag = feats.groupby("group")["knee_flag"].agg(["sum", "size"])
    print("knee_flag: " + ", ".join(f"{g} {int(r['sum'])}/{int(r['size'])}" for g, r in flag.loc[["TD", "AS"]].iterrows()))


if __name__ == "__main__":
    sys.exit(main())
