"""src/analysis.py: bootstrap reproducibility, median recovery, CI coverage on a synthetic set
with a known fixed - knee offset taken from the synthetic validation. No test reads data/."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from analysis import boot_ci, cluster_boot_ci, primary

ROOT = Path(__file__).resolve().parent.parent
CFG = {"analysis": {"bootstrap_n": 10000, "random_seed": 20261029, "ci_level": 0.95}}


def _frame(knee, fixed, group="TD"):
    n = len(knee)
    return pd.DataFrame({"label": [f"{group}{i}" for i in range(n)], "group": group, "id": np.arange(n),
                         "knee_flag": False, "periodic_delta_peak_knee": knee, "periodic_delta_peak_fixed": fixed})


def test_bootstrap_ci_identical_on_two_runs():
    x = np.random.default_rng(1).normal(0.3, 0.1, 185)
    assert boot_ci(x, np.median, 10000, 20261029, 0.95) == boot_ci(x, np.median, 10000, 20261029, 0.95)
    clusters = np.repeat(np.arange(62), 3)[:len(x)]
    assert (cluster_boot_ci(x, clusters, np.median, 2000, 20261029, 0.95)
            == cluster_boot_ci(x, clusters, np.median, 2000, 20261029, 0.95))
    # the seed matters (checked on the mean: bootstrap medians of 185 values sit on order statistics,
    # so two seeds can legitimately give identical median CIs)
    assert boot_ci(x, np.mean, 10000, 20261029, 0.95) != boot_ci(x, np.mean, 10000, 20261030, 0.95)


def test_median_of_constructed_paired_differences_is_recovered():
    d = np.array([-0.2, 0.05, 0.1, 0.25, 0.25, 0.3, 0.4, 0.9, 1.5])  # median 0.25
    knee = np.random.default_rng(2).uniform(0, 1, len(d))
    out = primary(_frame(knee, knee + d), CFG)
    row = out[(out["subset"] == "all") & (out["group"] == "TD")].iloc[0]
    assert row["median_d"] == pytest.approx(0.25, abs=1e-12)
    assert row["frac_d_gt_0"] == pytest.approx(8 / 9)
    assert row["mean_d"] == pytest.approx(d.mean())


def test_ci_covers_known_offset_from_synthetic_validation():
    sim = pd.read_csv(ROOT / "tables" / "sim_summary.csv")
    r = sim[(sim["condition"] == "a_baseline") & (sim["nlv"] == 0.02) & (sim["mode"] == "fixed-knee")].iloc[0]
    offset = float(r["periodic_delta_peak_median"])                       # 0.32072
    sd = float(r["periodic_delta_peak_q75"] - r["periodic_delta_peak_q25"]) / 1.349  # IQR -> normal SD
    rng = np.random.default_rng(20261029)
    knee = rng.uniform(0, 0.1, 185)
    noise = rng.normal(0, sd, 185)
    noise -= np.median(noise)  # symmetric noise with median exactly 0, so the true median d is the offset
    out = primary(_frame(knee, knee + offset + noise), CFG)
    row = out[(out["subset"] == "all") & (out["group"] == "TD")].iloc[0]
    assert row["median_ci_lower"] <= offset <= row["median_ci_upper"]
    assert row["median_d"] == pytest.approx(offset, abs=1e-12)
