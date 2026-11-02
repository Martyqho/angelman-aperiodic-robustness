"""Load config.yaml, with the types specparam needs."""

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = ROOT / "config.yaml"


def load_config(path=DEFAULT_CONFIG):
    """Return config.yaml as a dict.

    specparam.max_n_peaks is coerced with float(): YAML reads `.inf` as a float but
    a plain `inf` as the string 'inf', and float() accepts either.
    """
    with open(path) as f:
        cfg = yaml.safe_load(f)
    sp = cfg["specparam"]
    sp["max_n_peaks"] = float(sp["max_n_peaks"])
    sp["peak_threshold"] = float(sp["peak_threshold"])
    sp["min_peak_height"] = float(sp["min_peak_height"])
    sp["peak_width_limits_hz"] = [float(x) for x in sp["peak_width_limits_hz"]]
    sp["fit_range_hz"] = [float(x) for x in sp["fit_range_hz"]]
    return cfg
