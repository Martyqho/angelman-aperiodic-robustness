"""QC outputs for the fit pass (PROTOCOL gate items 2-5).

Reads data/features.csv and data/fits.npz (python src/fit_spectra.py first).

Writes
  tables/qc_fit_quality.csv  R^2 and fit-error quantiles per mode x group (aggregate)
  tables/qc_counts.csv       convergence and knee-flag counts by group, and by group x age bin
  data/qc/knee_flagged.csv   the knee-flagged recordings, listed (local)
  data/qc/hist_r2_error.png  R^2 and error histograms, per group, both modes (local)
  data/qc/overlay_selection.csv, data/qc/overlay_<label>.png, data/qc/overlays_<group>.png

Overlay selection (declared in NOTES.md 2026-11-04 before this ran), per group:
named case (TD133 / AS154), then the highest raw delta, then the knee-flagged recording
with the lowest knee-mode R^2 (fallback: knee frequency farthest from the group median,
log scale); a rule that hits an already-chosen recording takes the next in its ordering.
Then 7 drawn with RandomState(20261029), freshly seeded per group (TD first), from the
group's remaining labels in Fig 1 order.

Nothing here summarizes fixed-mode periodic delta or the fixed - knee difference.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import load_config

ROOT = Path(__file__).resolve().parent.parent
DATA, TABLES = ROOT / "data", ROOT / "tables"
QC = DATA / "qc"

GROUPS = ("TD", "AS")
MODES = ("knee", "fixed")
NAMED = {"TD": ("TD133", "footnoted offset"), "AS": ("AS154", "sharp 2.5 Hz peak")}
AGE_BINS = [0, 3, 6, 9, 12, 16]
QUANTILES = [0, 0.05, 0.25, 0.5, 0.75, 0.95, 1]
N_RANDOM = 7

# Same palette as the synthetic figure (dataviz reference palette, validated pair).
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
MODE_COLOR = {"knee": "#2a78d6", "fixed": "#eb6834"}
GROUP_ALPHA = {"TD": 0.45, "AS": 0.45}


def fit_quality(feats):
    rows = []
    for mode in MODES:
        for g in GROUPS:
            sub = feats[feats["group"] == g]
            for metric in ("r2", "error"):
                v = sub[f"{metric}_{mode}"].dropna()
                row = {"mode": mode, "group": g, "metric": metric, "n": len(v)}
                row.update({f"q{int(q * 100):02d}": v.quantile(q) for q in QUANTILES})
                row["mean"] = v.mean()
                rows.append(row)
    return pd.DataFrame(rows)


def counts(feats):
    f = feats.assign(age_bin=pd.cut(feats["age"], AGE_BINS, right=False).astype(str))
    rows = []
    for by in (["group"], ["group", "age_bin"]):
        for key, sub in f.groupby(by, sort=False):
            key = key if isinstance(key, tuple) else (key,)
            rows.append({"group": key[0], "age_bin": key[1] if len(key) > 1 else "all", "n": len(sub),
                         "converged_knee": int(sub["success_knee"].sum()),
                         "converged_fixed": int(sub["success_fixed"].sum()),
                         "knee_flag": int(sub["knee_flag"].sum()),
                         "knee_le_0": int((sub["knee_knee"] <= 0).sum()),
                         "knee_freq_out_of_range": int(((sub["knee_knee"] > 0) & sub["knee_flag"]).sum())})
    out = pd.DataFrame(rows)
    out["_g"] = out["group"].map({g: i for i, g in enumerate(GROUPS)})
    out["_b"] = out["age_bin"].map(lambda b: -1 if b == "all" else float(b.split(",")[0].strip("[")))
    return out.sort_values(["_g", "_b"]).drop(columns=["_g", "_b"]).reset_index(drop=True)


def select_overlays(feats, seed):
    chosen = []
    for g in GROUPS:
        sub = feats[feats["group"] == g]
        picked = []

        def take(ordered, reason):
            for lab in ordered:
                if lab not in picked:
                    picked.append(lab)
                    chosen.append({"group": g, "label": lab, "reason": reason})
                    return

        take([NAMED[g][0]], f"named: {NAMED[g][1]}")
        take(sub.sort_values("raw_delta", ascending=False)["label"], "highest raw delta")
        flagged = sub[sub["knee_flag"]]
        if len(flagged):
            take(flagged.sort_values("r2_knee")["label"], "knee-flagged, lowest knee R²")
        else:
            lk = np.log(sub["knee_freq_knee"])
            take(sub.assign(d=(lk - lk.median()).abs()).sort_values("d", ascending=False)["label"],
                 "knee frequency farthest from median")
        remaining = [lab for lab in sub["label"] if lab not in picked]
        for lab in np.random.RandomState(seed).choice(remaining, N_RANDOM, replace=False):
            picked.append(lab)
            chosen.append({"group": g, "label": lab, "reason": "random"})
    return pd.DataFrame(chosen)


def _style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=7, length=0)
    ax.grid(color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)


def _overlay_axes(ax_l, ax_r, fits, i, row, reason, legend=False):
    f = fits["freqs"]
    for ax in (ax_l, ax_r):
        _style(ax)
        ax.set_xscale("log")
        ax.set_xticks([1, 2, 4, 10, 20, 50], ["1", "2", "4", "10", "20", "50"])
        ax.axvspan(2, 4, color=MUTED, alpha=0.10, linewidth=0, zorder=0)
    ax_l.plot(f, fits["log_psd"][i], color=INK, linewidth=1.2, label="log10 PSD")
    for m in MODES:
        ax_l.plot(f, fits[f"aperiodic_{m}"][i], color=MODE_COLOR[m], linewidth=1.6, label=f"aperiodic, {m}")
        ax_r.plot(f, fits[f"residual_{m}"][i], color=MODE_COLOR[m], linewidth=0.9, alpha=0.55,
                  label=f"raw residual, {m}")
        ax_r.plot(f, fits[f"peak_{m}"][i], color=MODE_COLOR[m], linewidth=1.8, label=f"peak model, {m}")
    ax_r.axhline(0, color=AXIS, linewidth=0.8, zorder=1)
    kf = row["knee_freq_knee"]
    kf_txt = f"{kf:.2f} Hz" if np.isfinite(kf) else "none"
    ax_l.set_title(f"{row['label']}  age {row['age']:.1f}  ({reason})", loc="left", fontsize=8, color=INK)
    ax_r.set_title(f"R² knee {row['r2_knee']:.3f}, fixed {row['r2_fixed']:.3f}; knee freq {kf_txt}",
                   loc="left", fontsize=8, color=INK2)
    ax_l.set_ylabel("log10 power", fontsize=7, color=INK2)
    ax_r.set_ylabel("log10 ratio", fontsize=7, color=INK2)
    if legend:
        for ax in (ax_l, ax_r):
            ax.legend(frameon=False, fontsize=6.5, labelcolor=INK2, loc="upper right")


def plot_overlays(feats, fits, selection):
    idx = {lab: i for i, lab in enumerate(fits["labels"])}
    rows = feats.set_index("label")
    for _, s in selection.iterrows():
        fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.4), facecolor=SURFACE)
        _overlay_axes(a, b, fits, idx[s["label"]], rows.loc[s["label"]].to_dict() | {"label": s["label"]},
                      s["reason"], legend=True)
        a.set_xlabel("Hz (shaded 2–4)", fontsize=7, color=INK2)
        b.set_xlabel("Hz (shaded 2–4)", fontsize=7, color=INK2)
        fig.tight_layout()
        fig.savefig(QC / f"overlay_{s['label']}.png", dpi=130, facecolor=SURFACE)
        plt.close(fig)
    for g in GROUPS:
        sel = selection[selection["group"] == g].reset_index(drop=True)
        fig, axes = plt.subplots(len(sel), 2, figsize=(11, 2.3 * len(sel)), facecolor=SURFACE)
        for k, s in sel.iterrows():
            _overlay_axes(axes[k, 0], axes[k, 1], fits, idx[s["label"]],
                          rows.loc[s["label"]].to_dict() | {"label": s["label"]}, s["reason"], legend=(k == 0))
        fig.suptitle(f"{g}: raw PSD and aperiodic fits (left); raw residuals (thin) and peak models (thick), "
                     "knee blue, fixed orange (right)", x=0.01, ha="left", fontsize=9, color=INK2)
        fig.tight_layout(rect=(0, 0, 1, 0.985))
        fig.savefig(QC / f"overlays_{g}.png", dpi=110, facecolor=SURFACE)
        plt.close(fig)


def plot_hist(feats):
    fig, axes = plt.subplots(2, 2, figsize=(10, 6), facecolor=SURFACE)
    for c, g in enumerate(GROUPS):
        sub = feats[feats["group"] == g]
        for r, (metric, label) in enumerate([("r2", "R²"), ("error", "fit error (MAE, log10)")]):
            ax = axes[r, c]
            _style(ax)
            allv = feats[[f"{metric}_{m}" for m in MODES]].to_numpy().ravel()
            bins = np.linspace(np.nanmin(allv), np.nanmax(allv), 36)
            for m in MODES:
                ax.hist(sub[f"{metric}_{m}"].dropna(), bins=bins, color=MODE_COLOR[m], alpha=0.55,
                        label=f"{m} fit", edgecolor=SURFACE, linewidth=0.5)
            ax.set_title(f"{g}: {label}", loc="left", fontsize=9, color=INK)
            ax.set_ylabel("recordings", fontsize=8, color=INK2)
            if r == 0 and c == 0:
                ax.legend(frameon=False, fontsize=8, labelcolor=INK2)
    fig.tight_layout()
    fig.savefig(QC / "hist_r2_error.png", dpi=130, facecolor=SURFACE)
    plt.close(fig)


def main():
    cfg = load_config()
    seed = int(cfg["analysis"]["random_seed"])
    feats = pd.read_csv(DATA / "features.csv")
    with np.load(DATA / "fits.npz") as z:
        fits = {k: z[k] for k in z.files}
    if list(fits["labels"]) != list(feats["label"]):
        raise RuntimeError("features.csv and fits.npz disagree on order")
    QC.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(exist_ok=True)

    fit_quality(feats).to_csv(TABLES / "qc_fit_quality.csv", index=False, float_format="%.5g")
    counts(feats).to_csv(TABLES / "qc_counts.csv", index=False)
    feats.loc[feats["knee_flag"], ["label", "group", "age", "knee_knee", "knee_freq_knee", "exponent_knee",
                                   "r2_knee", "r2_fixed"]].to_csv(QC / "knee_flagged.csv", index=False,
                                                                  float_format="%.6g")
    plot_hist(feats)
    selection = select_overlays(feats, seed)
    selection.to_csv(QC / "overlay_selection.csv", index=False)
    plot_overlays(feats, fits, selection)


if __name__ == "__main__":
    main()
