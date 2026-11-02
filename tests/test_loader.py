"""One test per loader assertion (each one a trap found on Oct 25-26)."""

import logging
from pathlib import Path

import numpy as np
import pytest

import load_supplement as ls

XLSX = Path(__file__).resolve().parent.parent / "refs" / "supplementary-data-1.xlsx"


@pytest.fixture(scope="session")
def spectra():
    return ls.load_spectra(XLSX)


@pytest.fixture(scope="session")
def reported():
    return ls.load_reported(XLSX)


@pytest.fixture(scope="session")
def table():
    return ls.load_table(XLSX)


@pytest.fixture(scope="session")
def mapping(table, reported):
    return ls.map_table_to_spectra(table, reported)


# load_spectra

def test_spectra_98_freqs_1_to_49_5_step_0_5(spectra):
    freqs, psd, _ = spectra
    assert freqs.shape == (98,)
    assert freqs[0] == 1.0 and freqs[-1] == 49.5
    assert np.allclose(np.diff(freqs), 0.5)
    assert psd.shape == (344, 98)


def test_spectra_labels_td_then_as_skipping_separators(spectra):
    _, _, labels = spectra
    assert labels == [f"TD{i}" for i in range(1, 186)] + [f"AS{i}" for i in range(1, 160)]


def test_spectra_finite_and_positive(spectra):
    _, psd, _ = spectra
    assert np.isfinite(psd).all()
    assert (psd > 0).all()


# load_reported

def test_reported_strips_and_logs_star(caplog):
    with caplog.at_level(logging.WARNING, logger="load_supplement"):
        reported = ls.load_reported(XLSX)
    starred = [r.getMessage() for r in caplog.records if "*" in r.getMessage()]
    assert len(starred) == 1 and "TD133" in starred[0] and "offset" in starred[0]
    assert reported.loc["TD133", "offset"] == pytest.approx(-0.1314)


def test_reported_same_labels_same_order(reported, spectra):
    _, _, labels = spectra
    assert list(reported.index) == labels
    assert list(reported.columns) == [
        "delta", "theta", "alpha", "beta", "gamma", "total", "periodic_delta", "exponent", "offset",
    ]


# load_table

def test_table_344_rows(table):
    assert len(table) == 344


def test_table_185_nt_159_as(table):
    assert table["genotype"].value_counts().to_dict() == {"NT": 185, "AS": 159}


def test_table_95_unique_as_ids(table):
    assert table.loc[table["genotype"] == "AS", "id"].nunique() == 95


# map_table_to_spectra

def test_mapping_every_row_matched_within_tol(table, reported, mapping):
    matched = reported.loc[mapping.to_numpy(), ["exponent", "offset"]].to_numpy()
    dist = np.abs(table[["exponent", "offset"]].to_numpy() - matched).max(axis=1)
    assert mapping.notna().all()
    assert (dist <= 1e-5).all()
    with pytest.raises(ValueError, match="no Fig 1 match"):
        ls.map_table_to_spectra(table.assign(offset=table["offset"] + 1e-3), reported)


def test_mapping_344_unique_labels(mapping, reported):
    assert mapping.nunique() == 344
    assert set(mapping) == set(reported.index)


def test_mapping_group_agrees(table, mapping):
    assert (table["genotype"].map({"NT": "TD", "AS": "AS"}) == mapping.str[:2]).all()


def test_mapping_never_positional(table, reported, mapping):
    # Fig 2 is AS-first, Fig 1 TD-first ...
    assert mapping.iloc[0].startswith("AS") and reported.index[0].startswith("TD")
    # ... and a block swap is still wrong: within-block order differs.
    block_swap = [f"AS{i}" for i in range(1, 160)] + [f"TD{i}" for i in range(1, 186)]
    assert list(mapping) != block_swap
    # Row order of the table must not matter.
    shuffled = table.sample(frac=1, random_state=0)
    assert ls.map_table_to_spectra(shuffled, reported).sort_index().equals(mapping.sort_index())


# load_all

def test_load_all_tidy_aligned_with_psd(spectra, reported):
    tidy, psd = ls.load_all(XLSX)
    _, psd_ref, labels = spectra
    assert list(tidy["label"]) == labels
    assert np.array_equal(psd, psd_ref)
    assert list(tidy.columns[:6]) == ["label", "group", "id", "age", "sex", "del_status"]
    assert tidy["group"].value_counts().to_dict() == {"TD": 185, "AS": 159}
    assert tidy[["id", "age", "sex"]].notna().all().all()
    assert tidy.loc[tidy["group"] == "TD", "del_status"].isna().all()
    assert tidy.loc[tidy["group"] == "AS", "del_status"].isin(["del", "non-del", "unknown"]).all()
    assert np.array_equal(tidy[reported.columns].to_numpy(), reported.to_numpy())
    # TD1 is Fig 2 row 159 (the first NT row): age 4.25, Male.
    assert tidy.loc[0, ["age", "sex"]].tolist() == [4.25, "Male"]
