"""src/sim_gap_twins.py: the effects add up to M0 - T11, and the cells are built as declared."""

import numpy as np
import pytest

from sim_gap_twins import CELLS, effects, twin_cells


def test_effects_add_up_to_total():
    rng = np.random.default_rng(0)
    m = {c: rng.normal(0.2, 0.1, 50) for c in CELLS}
    e = effects(m)
    np.testing.assert_allclose(e["heterogeneity"] + e["peaks"] + e["residual"], e["total"], atol=1e-12)
    scalar = effects({c: float(v[0]) for c, v in m.items()})
    assert scalar["total"] == pytest.approx(m["M0"][0] - m["T11"][0])


def test_twin_cells_change_one_thing_at_a_time():
    rng = np.random.default_rng(1)
    ap, peak, resid, ap_med, noise = (rng.normal(size=98) for _ in range(5))
    c = twin_cells(ap, peak, resid, ap_med, noise)
    np.testing.assert_allclose(c["M0"] - c["T00"], ap_med - ap)
    np.testing.assert_allclose(c["T10"] - c["T00"], peak)
    np.testing.assert_allclose(c["T11"] - c["T01"], peak)
    np.testing.assert_allclose(c["T11"], ap + peak + resid)


def test_size_and_structure_add_up_to_residual_effect():
    rng = np.random.default_rng(2)
    e = effects({c: rng.normal(0.2, 0.1, 50) for c in CELLS})
    np.testing.assert_allclose(e["residual_size"] + e["residual_structure"], e["residual"], atol=1e-12)
    np.testing.assert_allclose((e["residual_size_no_peaks"] + e["residual_size_with_peaks"]) / 2,
                               e["residual_size"], atol=1e-12)


def test_size_matched_noise_has_the_residual_sd_and_t00_shape():
    rng = np.random.default_rng(3)
    ap, peak, ap_med = (rng.normal(size=98) for _ in range(3))
    resid, noise = rng.normal(0, 0.04, 98), rng.normal(0, 0.02, 98)
    c = twin_cells(ap, peak, resid, ap_med, noise)
    own = c["T0s"] - ap
    assert own.std() == pytest.approx(resid.std())
    np.testing.assert_allclose(own / own.std(), noise / noise.std())
    np.testing.assert_allclose(c["T1s"] - c["T0s"], peak)
