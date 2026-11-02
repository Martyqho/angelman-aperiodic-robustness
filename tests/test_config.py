import math

from config import DEFAULT_CONFIG, load_config


def test_max_n_peaks_is_float_inf():
    sp = load_config()["specparam"]
    assert isinstance(sp["max_n_peaks"], float) and math.isinf(sp["max_n_peaks"])


def test_max_n_peaks_string_inf_still_coerced(tmp_path):
    text = DEFAULT_CONFIG.read_text().replace("max_n_peaks: .inf", "max_n_peaks: inf")
    p = tmp_path / "config.yaml"
    p.write_text(text)
    assert math.isinf(load_config(p)["specparam"]["max_n_peaks"])
