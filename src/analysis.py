"""Primary analysis and secondaries (a)-(e) (PROTOCOL v2.3).

Inputs: data/features.csv (wide: one row per recording, `_knee` / `_fixed` columns per
mode, as written by src/fit_spectra.py) and config.yaml (seed, bootstrap_n, ci_level).
One exception, stated: secondary (a) and the gap decomposition in (d) compare the refit
with the paper's reported values, which are not in features.csv; those two read them
with load_supplement.load_reported. Nothing else is read; nothing is refit.

Primary. d = periodic delta (band rule over the fitted peak model, log10) under fixed
minus under knee, per recording. Median d overall and by group, recording-level
bootstrap of the median (bootstrap_n resamples, percentile CI at ci_level), fraction with
d > 0; mean d with bootstrap CI, labelled descriptive.
Sensitivities (predeclared): (1) without knee-flagged recordings; (2) cluster bootstrap
over participants (Fig 2 ID; key = group + id, all of a participant's recordings drawn
together). TD has 185 participants with one recording each, so its cluster CI equals its
recording CI by construction; the sensitivity bites in AS (95 participants) and overall.

Random numbers: every interval gets its own np.random.default_rng(seed), created fresh,
so each CI is reproducible on its own and does not depend on call order.

Labels, fixed before the first run (column `label` in every table):
  primary.csv                 predeclared (primary + both sensitivities); mean rows descriptive
  secondary_a.csv             predeclared, PROTOCOL (a); `note` carries the aggregation caveat
  secondary_b.csv             exploratory: delta-exponent and delta-offset (fixed - knee) by
                              group, as specified for this script. PROTOCOL (b) is AS-vs-TD
                              Cohen's d per mode; it is in group_contrast_by_mode.csv.
  secondary_c.csv             slopes on log2(age) (per doubling of age; PROTOCOL Statistics says
                              log age). TD per-mode slopes of exponent, offset and periodic delta
                              and their fixed - knee slope differences (= slope of the paired
                              difference): predeclared, PROTOCOL (c). Delta-exponent slopes with and
                              without knee-flagged recordings: predeclared (PROTOCOL expectations,
                              2026-11-05). AS rows: exploratory. CIs by paired bootstrap over
                              recordings. No GAM.
  secondary_d.csv             predeclared, PROTOCOL (d): knee-flag rate by group x age bin;
                              edge-rule count (peak-model delta == 0 and raw-residual delta > 0.3)
                              by mode and group; convergence, fit error and low-frequency peak
                              counts by mode and group. Gap decomposition (knee refit - reported,
                              src/gap_subsets.gap_subsets): post hoc, labelled so.
  secondary_e.csv             predeclared, PROTOCOL (e): the primary on the raw-residual delta.
  group_contrast_by_mode.csv  Welch Cohen's d (AS vs TD; d = mean difference / sqrt of the average
                              of the two variances) for periodic delta and exponent per mode with
                              bootstrap CI: predeclared, PROTOCOL (b). AS - TD median periodic
                              delta per mode with bootstrap CI, and the fixed-mode minus knee-mode
                              change in that contrast: exploratory (PROTOCOL names the contrast
                              only as Cohen's d).

Outputs (aggregates only): tables/primary.csv, tables/secondary_a.csv ... secondary_e.csv,
tables/group_contrast_by_mode.csv, figures/fig_primary_paired.png (pair-mean quantile bins,
not points), figures/fig_knee_flags_age.png, figures/fig_group_contrast_by_mode.png.
Local only: data/fig_primary_paired_points.png (per-recording scatter).
"""

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
GROUPS = ("TD", "AS")
MODES = ("knee", "fixed")
AGE_BINS = [0, 3, 6, 9, 12, 16]  # same bins as src/qc.py
EDGE_RESID_MIN = 0.3
AGGREGATION_NOTE = ("Fig 1e is not a single knee curve for any recording (src/fig1e_diagnostic.py); the reported "
                    "features likely come from fits on per-region (or per-channel) spectra, while the refit is one "
                    "fit of the averaged Fig 1a spectrum, so agreement is limited by aggregation level.")

# Dataviz reference palette: modes keep slots 1-2; groups take slots 3-4 (validated as a set
# of four; groups are below 3:1 contrast on the surface, so every group mark is direct-labelled).
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
MODE_COLOR = {"knee": "#2a78d6", "fixed": "#eb6834"}
GROUP_COLOR = {"TD": "#1baf7a", "AS": "#eda100"}


# ---------------------------------------------------------------- bootstrap

def _level(cfg):
    return float(cfg["analysis"]["ci_level"])


def _settings(cfg):
    return int(cfg["analysis"]["bootstrap_n"]), int(cfg["analysis"]["random_seed"]), _level(cfg)


def _pct(stats_, level):
    a = (1 - level) / 2
    return float(np.quantile(stats_, a)), float(np.quantile(stats_, 1 - a))


def boot_ci(x, stat, n_boot, seed, level):
    """Percentile CI of stat(x) over recording-level resamples."""
    x = np.asarray(x, float)
    idx = np.random.default_rng(seed).integers(0, len(x), (n_boot, len(x)))
    return _pct(stat(x[idx], axis=1), level)


def cluster_boot_ci(x, clusters, stat, n_boot, seed, level):
    """Percentile CI of stat(x) resampling clusters with replacement, all members together."""
    x = np.asarray(x, float)
    keys, inv = np.unique(np.asarray(clusters), return_inverse=True)
    members = [x[inv == k] for k in range(len(keys))]
    draws = np.random.default_rng(seed).integers(0, len(keys), (n_boot, len(keys)))
    out = np.array([stat(np.concatenate([members[j] for j in row])) for row in draws])
    return _pct(out, level)


def _clusters(df):
    return (df["group"].astype(str) + "_" + df["id"].astype(str)).to_numpy()


# ---------------------------------------------------------------- primary and (e)

def paired_summary(df, base, cfg):
    """Fixed - knee paired difference of column `base`, overall and per group, both subsets."""
    n_boot, seed, level = _settings(cfg)
    rows = []
    for subset, sub in (("all", df), ("no_knee_flag", df[~df["knee_flag"].astype(bool)])):
        for g in ("all",) + GROUPS:
            s = sub if g == "all" else sub[sub["group"] == g]
            if len(s) == 0:
                continue
            d = (s[f"{base}_fixed"] - s[f"{base}_knee"]).to_numpy(float)
            lo, hi = boot_ci(d, np.median, n_boot, seed, level)
            clo, chi = cluster_boot_ci(d, _clusters(s), np.median, n_boot, seed, level)
            mlo, mhi = boot_ci(d, np.mean, n_boot, seed, level)
            rows.append({"quantity": base, "subset": subset, "group": g, "n": len(d),
                         "n_participants": len(np.unique(_clusters(s))),
                         "median_d": float(np.median(d)), "median_ci_lower": lo, "median_ci_upper": hi,
                         "median_cluster_ci_lower": clo, "median_cluster_ci_upper": chi,
                         "frac_d_gt_0": float((d > 0).mean()),
                         "mean_d": float(d.mean()), "mean_ci_lower": mlo, "mean_ci_upper": mhi,
                         "median_knee": float(np.median(s[f"{base}_knee"])),
                         "median_fixed": float(np.median(s[f"{base}_fixed"])),
                         "label": ("predeclared: primary" if subset == "all" else
                                   "predeclared: sensitivity without knee-flagged")
                                  + "; cluster CI = predeclared participant sensitivity; mean = descriptive"})
    return pd.DataFrame(rows)


def primary(df, cfg):
    return paired_summary(df, "periodic_delta_peak", cfg)


def secondary_e(df, cfg):
    out = paired_summary(df, "periodic_delta_resid", cfg)
    out["label"] = out["label"].str.replace("predeclared: primary", "predeclared: PROTOCOL (e)")
    return out


# ---------------------------------------------------------------- (a)

def secondary_a(df, reported):
    from validation import agreement

    feats = {"exponent": ("exponent_knee", "exponent"), "offset": ("offset_knee", "offset"),
             "periodic_delta": ("periodic_delta_peak_knee", "periodic_delta")}
    j = df.set_index("label").join(reported[[r for _, r in feats.values()]].add_suffix("_reported"))
    rows = []
    for subset, sub in (("all", j), ("no_knee_flag", j[~j["knee_flag"].astype(bool)])):
        for g in GROUPS:
            s = sub[sub["group"] == g]
            for name, (ours, theirs) in feats.items():
                rows.append({"subset": subset, "group": g, "feature": name,
                             **agreement(s[ours].to_numpy(float), s[f"{theirs}_reported"].to_numpy(float)),
                             "label": "predeclared: PROTOCOL (a)", "note": AGGREGATION_NOTE})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- (b)

def secondary_b(df, cfg):
    n_boot, seed, level = _settings(cfg)
    rows = []
    flag = df["knee_flag"].astype(bool)
    for feat in ("exponent", "offset"):
        for subset, sub in (("all", df), ("knee_flagged", df[flag]), ("no_knee_flag", df[~flag])):
          for g in ("all",) + GROUPS:
            s = sub if g == "all" else sub[sub["group"] == g]
            if len(s) == 0:
                continue
            d = (s[f"{feat}_fixed"] - s[f"{feat}_knee"]).to_numpy(float)
            lo, hi = boot_ci(d, np.median, n_boot, seed, level)
            rows.append({"feature": f"delta_{feat}", "subset": subset, "group": g, "n": len(d),
                         "median_d": float(np.median(d)),
                         "median_ci_lower": lo, "median_ci_upper": hi, "mean_d": float(d.mean()),
                         "frac_d_gt_0": float((d > 0).mean()),
                         "median_knee": float(np.median(s[f"{feat}_knee"])),
                         "median_fixed": float(np.median(s[f"{feat}_fixed"])),
                         "label": ("exploratory (fixed - knee by group); PROTOCOL (b) is in group_contrast_by_mode.csv"
                                   if subset == "all" or feat == "offset" else
                                   "checks the 2026-11-05 expectation (delta-exponent ~ 0 when knee-flagged)")})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- (c)

def _slopes(x, y):
    """OLS slope of y on x along the last axis."""
    xc = x - x.mean(axis=-1, keepdims=True)
    return (xc * (y - y.mean(axis=-1, keepdims=True))).sum(axis=-1) / (xc ** 2).sum(axis=-1)


def secondary_c(df, cfg):
    n_boot, seed, level = _settings(cfg)
    rows = []
    for g in GROUPS:
        for subset, sub in (("all", df), ("no_knee_flag", df[~df["knee_flag"].astype(bool)])):
            s = sub[sub["group"] == g]
            x = np.log2(s["age"].to_numpy(float))
            idx = np.random.default_rng(seed).integers(0, len(x), (n_boot, len(x)))
            for feat in ("exponent", "offset", "periodic_delta_peak"):
                ys = {m: s[f"{feat}_{m}"].to_numpy(float) for m in MODES}
                ys["fixed_minus_knee"] = ys["fixed"] - ys["knee"]
                for what, y in ys.items():
                    if subset == "no_knee_flag" and not (feat == "exponent" and what == "fixed_minus_knee"):
                        continue  # without-flag rows only for delta-exponent (the declared expectation)
                    est = float(_slopes(x, y))
                    lo, hi = _pct(_slopes(x[idx], y[idx]), level)
                    if g == "AS":
                        label = "exploratory (AS)"
                    elif feat == "exponent" and what == "fixed_minus_knee":
                        label = "predeclared: PROTOCOL (c) slope difference; expectation 2026-11-05"
                    else:
                        label = "predeclared: PROTOCOL (c)"
                    rows.append({"group": g, "subset": subset, "feature": feat, "series": what, "n": len(x),
                                 "slope_per_doubling_of_age": est, "ci_lower": lo, "ci_upper": hi,
                                 "label": label})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- (d)

def secondary_d(df, reported):
    from gap_subsets import gap_subsets

    rows = []

    def add(block, group, mode, stratum, measure, value, n, label):
        rows.append({"block": block, "group": group, "mode": mode, "stratum": stratum, "measure": measure,
                     "value": value, "n": n, "label": label})

    f = df.assign(age_bin=pd.cut(df["age"], AGE_BINS, right=False).astype(str))
    for g in GROUPS:
        s = f[f["group"] == g]
        for b in ["all"] + [str(c) for c in pd.cut(pd.Series([0.0]), AGE_BINS, right=False).cat.categories]:
            sb = s if b == "all" else s[s["age_bin"] == b]
            add("knee_flag_rate", g, "knee", b, "flagged", int(sb["knee_flag"].sum()), len(sb), "predeclared: PROTOCOL (d)")
            add("knee_flag_rate", g, "knee", b, "rate", float(sb["knee_flag"].mean()), len(sb), "predeclared: PROTOCOL (d)")
        for m in MODES:
            edge = (s[f"periodic_delta_peak_{m}"] == 0) & (s[f"periodic_delta_resid_{m}"] > EDGE_RESID_MIN)
            add("edge_rule_count", g, m, "all", "peak_delta_0_and_resid_delta_gt_0.3", int(edge.sum()), len(s),
                "predeclared: PROTOCOL (d)")
            add("diagnostics", g, m, "all", "converged", int(s[f"success_{m}"].sum()), len(s), "predeclared: PROTOCOL (d)")
            add("diagnostics", g, m, "all", "median_fit_error", float(s[f"error_{m}"].median()), len(s),
                "predeclared: PROTOCOL (d)")
            add("diagnostics", g, m, "all", "median_n_peaks", float(s[f"n_peaks_{m}"].median()), len(s),
                "predeclared: PROTOCOL (d)")
            add("diagnostics", g, m, "all", "median_n_peaks_below_6hz", float(s[f"n_peaks_below_6hz_{m}"].median()),
                len(s), "predeclared: PROTOCOL (d)")
            add("diagnostics", g, m, "all", "mean_n_peaks_below_6hz", float(s[f"n_peaks_below_6hz_{m}"].mean()),
                len(s), "predeclared: PROTOCOL (d)")

    j = df.set_index("label")
    gd = pd.DataFrame({"group": j["group"], "pd_refit": j["periodic_delta_peak_knee"],
                       "pd_reported": reported.loc[j.index, "periodic_delta"]})
    for _, r in gap_subsets(gd).iterrows():
        for measure in ("mean_gap", "median_gap", "share_of_group_gap"):
            add("gap_decomposition", r["group"], "knee", r["subset"], measure, r[measure], int(r["n"]),
                "post hoc (2026-11-14): knee refit - reported periodic delta")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- group contrast

def _cohen_d(a, b, axis=-1):
    return (a.mean(axis=axis) - b.mean(axis=axis)) / np.sqrt((a.var(axis=axis, ddof=1) + b.var(axis=axis, ddof=1)) / 2)


def group_contrast(df, cfg):
    n_boot, seed, level = _settings(cfg)
    as_, td = df[df["group"] == "AS"], df[df["group"] == "TD"]
    rows = []
    for feat in ("periodic_delta_peak", "exponent"):
        rng = np.random.default_rng(seed)
        ia = rng.integers(0, len(as_), (n_boot, len(as_)))
        it = rng.integers(0, len(td), (n_boot, len(td)))
        med = {}
        for m in MODES:
            a, t = as_[f"{feat}_{m}"].to_numpy(float), td[f"{feat}_{m}"].to_numpy(float)
            lo, hi = _pct(_cohen_d(a[ia], t[it]), level)
            rows.append({"feature": feat, "mode": m, "measure": "welch_cohens_d_AS_vs_TD",
                         "estimate": float(_cohen_d(a, t)), "ci_lower": lo, "ci_upper": hi,
                         "median_AS": float(np.median(a)), "median_TD": float(np.median(t)),
                         "n_AS": len(a), "n_TD": len(t), "label": "predeclared: PROTOCOL (b)"})
            if feat == "periodic_delta_peak":
                med[m] = np.median(a[ia], axis=1) - np.median(t[it], axis=1)
                lo, hi = _pct(med[m], level)
                rows.append({"feature": feat, "mode": m, "measure": "median_AS_minus_median_TD",
                             "estimate": float(np.median(a) - np.median(t)), "ci_lower": lo, "ci_upper": hi,
                             "median_AS": float(np.median(a)), "median_TD": float(np.median(t)),
                             "n_AS": len(a), "n_TD": len(t), "label": "exploratory"})
        if feat == "periodic_delta_peak":
            est = {m: float(np.median(as_[f"{feat}_{m}"]) - np.median(td[f"{feat}_{m}"])) for m in MODES}
            lo, hi = _pct(med["fixed"] - med["knee"], level)
            rows.append({"feature": feat, "mode": "fixed_minus_knee",
                         "measure": "change_in_median_AS_minus_median_TD", "estimate": est["fixed"] - est["knee"],
                         "ci_lower": lo, "ci_upper": hi, "median_AS": np.nan, "median_TD": np.nan,
                         "n_AS": len(as_), "n_TD": len(td), "label": "exploratory (paired over the same resamples)"})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- figures

def _style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8, length=0)
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def fig_primary_paired(df, prim, out, n_bins=5):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), facecolor=SURFACE, sharey=True)
    for ax, g in zip(axes, GROUPS):
        _style(ax)
        s = df[df["group"] == g]
        d = (s["periodic_delta_peak_fixed"] - s["periodic_delta_peak_knee"]).to_numpy(float)
        pm = ((s["periodic_delta_peak_fixed"] + s["periodic_delta_peak_knee"]) / 2).to_numpy(float)
        r = prim[(prim["subset"] == "all") & (prim["group"] == g)].iloc[0]
        ax.axhspan(r["median_ci_lower"], r["median_ci_upper"], color=GROUP_COLOR[g], alpha=0.18, linewidth=0)
        ax.axhline(r["median_d"], color=GROUP_COLOR[g], linewidth=2)
        ax.axhline(0, color=AXIS, linewidth=1)
        bins = pd.qcut(pm, n_bins, labels=False, duplicates="drop")
        for b in np.unique(bins):
            m = bins == b
            q25, q50, q75 = np.percentile(d[m], [25, 50, 75])
            x = np.median(pm[m])
            ax.plot([x, x], [q25, q75], color=INK2, linewidth=2, solid_capstyle="round")
            ax.plot(x, q50, "o", ms=8, color=GROUP_COLOR[g], mec=SURFACE, mew=2)
            ax.text(x, q75, f"n={m.sum()}", ha="center", va="bottom", fontsize=7, color=MUTED)
        ax.text(0.02, 0.97, f"{g} median d {r['median_d']:+.3f} (line)\n95% CI {r['median_ci_lower']:+.3f} to "
                f"{r['median_ci_upper']:+.3f} (band)", transform=ax.transAxes, ha="left", va="top",
                fontsize=8, color=INK)
        ax.set_title(f"{g} (n = {len(d)})", loc="left", fontsize=10, color=INK)
        ax.set_xlabel("pair mean of periodic delta, (fixed + knee) / 2, quintile bins", fontsize=8, color=INK2)
    axes[0].set_ylabel("d = fixed − knee periodic delta (log10)", fontsize=8, color=INK2)
    fig.suptitle("Primary: paired fixed − knee periodic delta. Dots: bin median; bars: bin IQR; line and band: "
                 "group median and recording-bootstrap 95% CI.", x=0.01, ha="left", fontsize=9, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_primary_points(df, out):
    fig, ax = plt.subplots(figsize=(7, 4.5), facecolor=SURFACE)
    _style(ax)
    for g in GROUPS:
        s = df[df["group"] == g]
        d = s["periodic_delta_peak_fixed"] - s["periodic_delta_peak_knee"]
        pm = (s["periodic_delta_peak_fixed"] + s["periodic_delta_peak_knee"]) / 2
        ax.scatter(pm, d, s=14, color=GROUP_COLOR[g], alpha=0.7, edgecolor=SURFACE, linewidth=0.5, label=g)
    ax.axhline(0, color=AXIS, linewidth=1)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2)
    ax.set_xlabel("pair mean", fontsize=8, color=INK2)
    ax.set_ylabel("fixed − knee", fontsize=8, color=INK2)
    ax.set_title("Local only: per-recording paired differences (not for publication)", loc="left", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=130, facecolor=SURFACE)
    plt.close(fig)


def fig_knee_flags_age(sec_d, out):
    k = sec_d[(sec_d["block"] == "knee_flag_rate") & (sec_d["stratum"] != "all")]
    bins = list(dict.fromkeys(k["stratum"]))
    fig, ax = plt.subplots(figsize=(8, 3.8), facecolor=SURFACE)
    _style(ax)
    off = {"TD": -0.12, "AS": 0.12}
    for g in GROUPS:
        r = k[(k["group"] == g) & (k["measure"] == "rate")].set_index("stratum").loc[bins]
        c = k[(k["group"] == g) & (k["measure"] == "flagged")].set_index("stratum").loc[bins]
        xs = np.arange(len(bins)) + off[g]
        ax.vlines(xs, 0, r["value"], color=GROUP_COLOR[g], linewidth=3)
        ax.plot(xs, r["value"], "o", ms=8, color=GROUP_COLOR[g], mec=SURFACE, mew=2, label=g)
        for x, v, f_, n in zip(xs, r["value"], c["value"], c["n"]):
            ax.text(x, v + 0.012, f"{g}\n{int(f_)}/{int(n)}", ha="center", va="bottom", fontsize=7, color=INK2)
    ax.set_xticks(range(len(bins)), [b.replace(", ", "–").strip("[)") + " y" for b in bins])
    ax.set_ylabel("knee-flagged fraction", fontsize=8, color=INK2)
    ax.set_ylim(0, max(0.35, k[k["measure"] == "rate"]["value"].max() + 0.08))
    ax.legend(frameon=False, fontsize=8, labelcolor=INK2, loc="upper right")
    ax.set_title("Knee flags (knee ≤ 0 or knee frequency outside 0.5–20 Hz) by age bin", loc="left", fontsize=10,
                 color=INK)
    fig.tight_layout()
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def fig_group_contrast(gc, out):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.8), facecolor=SURFACE)
    for ax in (a1, a2):
        _style(ax)
        ax.axhline(0, color=AXIS, linewidth=1)
    med = gc[gc["measure"] == "median_AS_minus_median_TD"].set_index("mode")
    for i, m in enumerate(MODES):
        r = med.loc[m]
        a1.plot([i, i], [r["ci_lower"], r["ci_upper"]], color=MODE_COLOR[m], linewidth=2.5, solid_capstyle="round")
        a1.plot(i, r["estimate"], "o", ms=9, color=MODE_COLOR[m], mec=SURFACE, mew=2)
        a1.text(i + 0.08, r["estimate"], f"{r['estimate']:+.3f}\n[{r['ci_lower']:+.3f}, {r['ci_upper']:+.3f}]",
                fontsize=8, color=INK, va="center")
    a1.set_xticks(range(len(MODES)), [f"{m} fit" for m in MODES])
    a1.set_xlim(-0.5, 1.9)
    a1.set_ylabel("AS − TD median periodic delta (log10)", fontsize=8, color=INK2)
    a1.set_title("AS − TD median difference (exploratory)", loc="left", fontsize=10, color=INK)
    cd = gc[gc["measure"] == "welch_cohens_d_AS_vs_TD"]
    for j, feat in enumerate(("periodic_delta_peak", "exponent")):
        for i, m in enumerate(MODES):
            r = cd[(cd["feature"] == feat) & (cd["mode"] == m)].iloc[0]
            x = j + (i - 0.5) * 0.3
            a2.plot([x, x], [r["ci_lower"], r["ci_upper"]], color=MODE_COLOR[m], linewidth=2.5, solid_capstyle="round")
            a2.plot(x, r["estimate"], "o", ms=9, color=MODE_COLOR[m], mec=SURFACE, mew=2,
                    label=f"{m} fit" if j == 0 else None)
            a2.text(x, r["ci_upper"], f"{r['estimate']:.2f}", ha="center", va="bottom", fontsize=8, color=INK)
    a2.set_xticks([0, 1], ["periodic delta", "exponent"])
    a2.set_xlim(-0.6, 1.6)
    a2.set_ylabel("Welch Cohen's d, AS vs TD", fontsize=8, color=INK2)
    a2.set_title("AS vs TD effect size (PROTOCOL (b))", loc="left", fontsize=10, color=INK)
    a2.legend(frameon=False, fontsize=8, labelcolor=INK2, loc="upper right")
    fig.suptitle("Group contrast under each fit mode. Dots: estimate; bars: bootstrap 95% CI.", x=0.01, ha="left",
                 fontsize=9, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)


# ---------------------------------------------------------------- main

def main():
    import load_supplement as ls
    from config import load_config

    logging.getLogger("load_supplement").setLevel(logging.ERROR)
    cfg = load_config()
    df = pd.read_csv(ROOT / "data" / "features.csv")
    if not (df["success_knee"].all() and df["success_fixed"].all()):
        df = df[df["success_knee"] & df["success_fixed"]]  # PROTOCOL Inclusion: both fits converge
    df["knee_flag"] = df["knee_flag"].astype(bool)
    reported = ls.load_reported(ROOT / cfg["paper"]["supplement_file"])  # (a) and (d) gap decomposition only

    tables, figures = ROOT / "tables", ROOT / "figures"
    tables.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)
    prim = primary(df, cfg)
    sec_d = secondary_d(df, reported)
    gc = group_contrast(df, cfg)
    outputs = {"primary": prim, "secondary_a": secondary_a(df, reported), "secondary_b": secondary_b(df, cfg),
               "secondary_c": secondary_c(df, cfg), "secondary_d": sec_d, "secondary_e": secondary_e(df, cfg),
               "group_contrast_by_mode": gc}
    for name, t in outputs.items():
        t.to_csv(tables / f"{name}.csv", index=False, float_format="%.5g")

    fig_primary_paired(df, prim, figures / "fig_primary_paired.png")
    fig_knee_flags_age(sec_d, figures / "fig_knee_flags_age.png")
    fig_group_contrast(gc, figures / "fig_group_contrast_by_mode.png")
    fig_primary_points(df, ROOT / "data" / "fig_primary_paired_points.png")


if __name__ == "__main__":
    main()
