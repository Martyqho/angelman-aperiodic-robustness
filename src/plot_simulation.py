"""figures/fig_synthetic_validation.png: the synthetic grid, both fit modes, both noise levels.

Rows: noise level. Columns: (1) raw delta by condition; (2) peak-model periodic delta
by condition against the baseline's 95% spread; (3) recovered vs injected periodic delta.
Marks: median dot with IQR bar; knee and fixed side by side.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reference palette (dataviz skill, references/palette.md), light mode; validated as a pair.
SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
MODE_COLOR = {"knee": "#2a78d6", "fixed": "#eb6834"}
MODE_LABEL = {"knee": "knee fit", "fixed": "fixed fit"}

CONDITIONS = [("a_baseline", 0.0, "(a)\nbaseline"), ("b_offset", 0.0, "(b)\noffset"),
              ("c_exponent", 0.0, "(c)\nexponent"), ("c_exponent_kf", 0.0, "(c′)\nexp., same\nknee freq"),
              ("d_peak", 1.0, "(d)\npeak"),
              ("e_all", 1.0, "(e)\nall three")]
OFFSET = {"knee": -0.13, "fixed": 0.13}


def _style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8, length=0)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def _dot_iqr(ax, x, values, mode):
    v = values[np.isfinite(values)]
    q25, med, q75 = np.percentile(v, [25, 50, 75])
    c = MODE_COLOR[mode]
    ax.plot([x, x], [q25, q75], color=c, linewidth=2, solid_capstyle="round", zorder=2)
    ax.plot(x, med, "o", ms=6, color=c, mec=SURFACE, mew=1.5, zorder=3)


def _by_condition(ax, sims, nlv, metric):
    for i, (cond, a_scale, _) in enumerate(CONDITIONS):
        g = sims[(sims["condition"] == cond) & (sims["a_scale"] == a_scale) & (sims["nlv"] == nlv)]
        for mode in MODE_COLOR:
            _dot_iqr(ax, i + OFFSET[mode], g.loc[g["mode"] == mode, metric].to_numpy(), mode)
    ax.set_xticks(range(len(CONDITIONS)), [c[2] for c in CONDITIONS])
    ax.set_xlim(-0.6, len(CONDITIONS) - 0.4)


def plot(sims, base, shifts, out):
    nlvs = sorted(sims["nlv"].unique())
    fig, axes = plt.subplots(len(nlvs), 3, figsize=(14, 4.4 * len(nlvs)), facecolor=SURFACE,
                             gridspec_kw={"width_ratios": [1.15, 1.15, 1]})
    axes = np.atleast_2d(axes)

    for r, nlv in enumerate(nlvs):
        ax1, ax2, ax3 = axes[r]
        for ax in axes[r]:
            _style(ax)

        # (1) raw delta: rises with offset, falls with exponent (pivot at 1 Hz), rises with a peak.
        _by_condition(ax1, sims, nlv, "raw_delta")
        ax1.set_yscale("log")
        ax1.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([10, 20, 50, 100, 200, 500]))
        ax1.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
        ax1.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        truth = (sims[(sims["condition"] == "a_baseline") & (sims["nlv"] == nlv)]["true_raw_delta"].iloc[0])
        ax1.axhline(truth, color=MUTED, linewidth=1, linestyle=(0, (3, 3)), zorder=1)
        ax1.set_ylabel("raw delta, 2–4 Hz (µV², log scale)", color=INK2, fontsize=9)

        # (2) peak-model periodic delta vs the no-peak baseline's 95% spread.
        a = sims[(sims["condition"] == "a_baseline") & (sims["nlv"] == nlv)]
        for mode, c in MODE_COLOR.items():
            v = a.loc[a["mode"] == mode, "periodic_delta_peak"]
            ax2.axhspan(v.quantile(0.025), v.quantile(0.975), color=c, alpha=0.10, linewidth=0, zorder=0)
        _by_condition(ax2, sims, nlv, "periodic_delta_peak")
        ax2.axhline(shifts["d_periodic_delta"], color=MUTED, linewidth=1, linestyle=(0, (3, 3)), zorder=1)
        ax2.text(len(CONDITIONS) - 0.45, shifts["d_periodic_delta"], "injected", color=INK2, fontsize=8,
                 ha="right", va="bottom")
        ax2.set_ylabel("periodic delta, peak model (log10)", color=INK2, fontsize=9)

        # (3) recovery: recovered minus the baseline median, vs injected, peak-only and all-three.
        lim = 0
        for cond, marker, name in [("d_peak", "o", "peak only"), ("e_all", "s", "all three")]:
            g = sims[(sims["condition"] == cond) & (sims["nlv"] == nlv)]
            for mode, c in MODE_COLOR.items():
                base_med = a.loc[a["mode"] == mode, "periodic_delta_peak"].median()
                xs, meds, lo, hi = [], [], [], []
                for _, h in g[g["mode"] == mode].groupby("a_scale"):
                    v = h["periodic_delta_peak"].to_numpy() - base_med
                    q25, m, q75 = np.nanpercentile(v, [25, 50, 75])
                    xs.append(h["true_periodic_delta"].iloc[0]); meds.append(m); lo.append(m - q25); hi.append(q75 - m)
                jitter = (0.012 if mode == "fixed" else -0.012) + (0.006 if cond == "e_all" else -0.006)
                xs = np.array(xs) + jitter
                ax3.errorbar(xs, meds, yerr=[lo, hi], fmt=marker, ms=6, color=c, mec=SURFACE, mew=1.2,
                             elinewidth=2, capsize=0, zorder=3,
                             label=f"{MODE_LABEL[mode]}, {name}" if r == 0 else None)
                lim = max(lim, max(xs), max(np.array(meds) + np.array(hi)))
        lim *= 1.08
        ax3.plot([0, lim], [0, lim], color=MUTED, linewidth=1, linestyle=(0, (3, 3)), zorder=1)
        ax3.text(lim * 0.1, lim * 0.1 - lim * 0.02, "identity", color=INK2, fontsize=8, ha="left", va="top")
        ax3.set_xlim(0, lim); ax3.set_ylim(min(0, ax3.get_ylim()[0]), lim)
        ax3.set_xlabel("injected periodic delta (log10)", color=INK2, fontsize=9)
        ax3.set_ylabel("recovered − baseline median (log10)", color=INK2, fontsize=9)

        ax1.set_title(f"Raw delta    noise {nlv:g}", loc="left", color=INK, fontsize=10)
        ax2.set_title(f"Periodic delta (paper's quantity)    noise {nlv:g}", loc="left", color=INK, fontsize=10)
        ax3.set_title(f"Recovery of the 3 Hz peak    noise {nlv:g}", loc="left", color=INK, fontsize=10)

    # One legend for the modes, one for the recovery markers; both on the first row.
    handles = [plt.Line2D([], [], marker="o", color=c, linewidth=2, ms=6, label=MODE_LABEL[m])
               for m, c in MODE_COLOR.items()]
    axes[0, 0].legend(handles=handles, frameon=False, fontsize=8, labelcolor=INK2, loc="upper left")
    axes[0, 2].legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc="upper left")

    fig.suptitle(
        f"Synthetic validation, {sims['replicate'].nunique()} seeds per cell. Baseline from TD knee fits: "
        f"b₀ = {base['b0']:.2f}, k₀ = {base['k0']:.1f}, χ₀ = {base['chi0']:.2f}. "
        f"Shifts from reported AS−TD medians: Δb = {shifts['db']:+.2f}, Δχ = {shifts['dchi']:+.2f}, "
        f"peak A = {shifts['A']:.2f} at 3 Hz, σ = 0.5.\n"
        "Dots: median; bars: IQR. Shaded: no-peak baseline 95% spread (per mode). Dashed: truth.",
        x=0.01, ha="left", fontsize=9, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)
