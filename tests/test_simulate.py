"""Generator and fitter sanity checks; if these fail, nothing downstream means anything."""

import numpy as np
import pytest
from specparam.sim import sim_power_spectrum

from config import load_config
from fit_spectra import fit_one
from simulate import FREQS, conditions, generate

# Baseline from the 185 TD knee fits (simulate.baseline_from_td, 2026-11-03).
B0, K0, CHI0 = 2.359707, 16.518708, 2.208844
A = 0.6031


@pytest.fixture(scope="module")
def cfg():
    return load_config()


def test_noise_free_baseline_recovers_b0_chi0(cfg):
    psd = generate(B0, K0, CHI0, A=0.0, nlv=0.0)
    res = fit_one(FREQS, psd, "knee", cfg)
    assert res["success"]
    assert res["exponent"] == pytest.approx(CHI0, abs=1e-3), res
    assert res["offset"] == pytest.approx(B0, abs=1e-3), res
    assert res["n_peaks"] == 0, res


@pytest.mark.parametrize("nlv", [0.0, 0.05])
def test_generator_matches_specparam_sim(nlv):
    seed = 20261029
    np.random.seed(seed)
    f, p = sim_power_spectrum([1, 49.5], [B0, K0, CHI0], [3.0, A, 0.5], nlv=nlv, freq_res=0.5)
    assert np.allclose(f, FREQS)
    assert np.allclose(np.log10(p), np.log10(generate(B0, K0, CHI0, A=A, nlv=nlv, seed=seed)),
                       atol=1e-12, rtol=0)


def test_c_prime_keeps_knee_frequency():
    base = {"b0": B0, "k0": K0, "chi0": CHI0}
    cells = {c: p for c, _, p in conditions(base, {"db": 0.842, "dchi": 0.393, "A": A})}
    kf0 = K0 ** (1 / CHI0)
    c, cp = cells["c_exponent"], cells["c_exponent_kf"]
    assert c["chi"] == cp["chi"] == pytest.approx(CHI0 + 0.393)
    assert cp["k"] ** (1 / cp["chi"]) == pytest.approx(kf0, rel=1e-12)
    assert c["k"] ** (1 / c["chi"]) < kf0 - 0.5  # (c) moves the knee down; (c') does not
