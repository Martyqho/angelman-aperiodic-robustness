"""Load Supplementary Data 1 of Montgomery et al. 2026 (Commun Med 6:532).

Four loaders, each of which checks the traps found on Oct 25-26 and raises
ValueError if the workbook does not look the way we expect:

  load_spectra(path)          'Fig 1a'        -> freqs (98,), psd (344, 98), labels
  load_reported(path)         'Fig 1b,d,f,g'  -> reported features, indexed by label
  load_table(path)            'Fig 2'         -> per-recording demographics table
  map_table_to_spectra(...)   Fig 2 row -> Fig 1 label, by nearest (exponent, offset)

load_all(path) joins them into one tidy table and returns it with the psd array.

Fig 2's Raw_delta and Periodic_delta columns are read but not used: they are on an
undocumented convention (see PROTOCOL.md, Conventions).
"""

import logging
import re
from functools import lru_cache
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

log = logging.getLogger(__name__)

N_TD, N_AS = 185, 159
N_REC = N_TD + N_AS
N_AS_SUBJECTS = 95
FREQS_EXPECTED = np.arange(1.0, 49.5 + 0.25, 0.5)  # 1.0, 1.5, ..., 49.5 Hz (98 bins)
LABELS_EXPECTED = [f"TD{i}" for i in range(1, N_TD + 1)] + [f"AS{i}" for i in range(1, N_AS + 1)]
LABEL_RE = re.compile(r"^(TD|AS)\d+$")

# 'Fig 1b,d,f,g' row label (column 1, normalized) -> feature name.
# Matched on the name only: the sheet's band ranges (e.g. 'Beta 12-30 Hz') and
# whitespace ('Aperiodic     offset', non-breaking spaces) are not reliable keys.
REPORTED_ROWS = {
    "delta": "delta",
    "theta": "theta",
    "alpha": "alpha",
    "beta": "beta",
    "gamma": "gamma",
    "total": "total",
    "periodic delta": "periodic_delta",
    "aperiodic exponent": "exponent",
    "aperiodic offset": "offset",
}
REPORTED_FEATURES = list(REPORTED_ROWS.values())

# 'Fig 2' header -> column name in load_table().
TABLE_COLUMNS = {
    "ID": "id",
    "Age": "age",
    "Sex": "sex",
    "Genotype": "genotype",  # 'NT' / 'AS' (Fig 1 calls NT 'TD')
    "del vs non-del": "del_status",  # 'del' / 'non-del' / 'unknown'; None for NT
    "Raw_delta": "raw_delta",  # unused: undocumented convention
    "log_10(Raw_delta)": "log10_raw_delta",  # unused
    "Periodic_delta": "periodic_delta",  # unused: undocumented convention
    "Aperiodic_exponent": "exponent",
    "Aperiodic_offset": "offset",
}
GENOTYPE_TO_GROUP = {"NT": "TD", "AS": "AS"}


@lru_cache(maxsize=None)
def _sheet_rows(path, sheet):
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    try:
        return tuple(wb[sheet].iter_rows(values_only=True))
    finally:
        wb.close()


def _rows(path, sheet):
    return _sheet_rows(str(Path(path).resolve()), sheet)


def _label_columns(header, start):
    """Return [(col, label)] for subject columns, skipping unlabeled separator columns."""
    cols = []
    for j, h in enumerate(header[start:], start=start):
        if h is None or str(h).strip() == "":
            continue
        h = str(h).strip()
        if not LABEL_RE.match(h):
            raise ValueError(f"unexpected column header {h!r} at column {j}")
        cols.append((j, h))
    return cols


def _check_labels(labels, what):
    if list(labels) != LABELS_EXPECTED:
        raise ValueError(f"{what}: expected labels TD1..TD{N_TD} then AS1..AS{N_AS}, got {len(labels)} labels")


def _norm(s):
    return " ".join(str(s).split("\n")[0].split()).casefold()


def load_freq_sheet(path, sheet):
    """Return (freqs (98,), values (344, 98), labels) from a Fig 1a/1c/1e-style sheet.

    values[i] is the row for labels[i]; all values are checked finite.
    """
    rows = _rows(path, sheet)
    h = next(i for i, r in enumerate(rows) if r and r[0] == "Frequency")
    cols = _label_columns(rows[h], start=1)
    labels = [lab for _, lab in cols]
    _check_labels(labels, sheet)

    body = [r for r in rows[h + 1:] if r and r[0] is not None]
    freqs = np.array([float(r[0]) for r in body])
    if freqs.shape != FREQS_EXPECTED.shape or not np.allclose(freqs, FREQS_EXPECTED):
        raise ValueError(f"{sheet}: expected 98 frequencies 1.0-49.5 Hz step 0.5, got {freqs.size}")

    values = np.array([[np.nan if r[j] is None else float(r[j]) for r in body] for j, _ in cols])
    if not np.isfinite(values).all():
        raise ValueError(f"{sheet}: non-finite values")
    return freqs, values, labels


def load_spectra(path):
    """Return (freqs (98,), psd (344, 98), labels) from sheet 'Fig 1a'.

    psd[i] is the spectrum of labels[i], linear uV^2/Hz.
    """
    freqs, psd, labels = load_freq_sheet(path, "Fig 1a")
    if not (psd > 0).all():
        raise ValueError("Fig 1a: non-positive PSD values")
    return freqs, psd, labels


def load_reported(path):
    """Return the paper's reported per-recording features from 'Fig 1b,d,f,g'.

    DataFrame indexed by label (TD1..TD185, AS1..AS159), columns REPORTED_FEATURES.
    A trailing '*' on a cell (e.g. TD133 offset '-0.1314*') is stripped and logged.
    """
    rows = _rows(path, "Fig 1b,d,f,g")
    h = next(i for i, r in enumerate(rows) if r and "TD1" in r)
    cols = _label_columns(rows[h], start=2)
    labels = [lab for _, lab in cols]
    _check_labels(labels, "Fig 1b,d,f,g")

    data = {}
    for r in rows[h + 1:]:
        if not r or r[1] is None:
            continue
        name = REPORTED_ROWS.get(_norm(r[1]))
        if name is None:
            continue
        if name in data:
            raise ValueError(f"Fig 1b,d,f,g: duplicate row for {name!r}")
        values = []
        for j, lab in cols:
            v = r[j]
            if isinstance(v, str):
                s = v.strip()
                if s.endswith("*"):
                    log.warning("Fig 1b,d,f,g: stripped trailing '*' from %s %s (%r)", lab, name, v)
                    s = s.rstrip("*").strip()
                v = s
            values.append(float(v))
        data[name] = values
    missing = [f for f in REPORTED_FEATURES if f not in data]
    if missing:
        raise ValueError(f"Fig 1b,d,f,g: missing rows {missing}")

    df = pd.DataFrame(data, index=pd.Index(labels, name="label"))[REPORTED_FEATURES]
    if not np.isfinite(df.to_numpy()).all():
        raise ValueError("Fig 1b,d,f,g: non-finite reported values")
    return df


def load_table(path):
    """Return sheet 'Fig 2' as a DataFrame, one row per recording (AS block first).

    Columns are renamed per TABLE_COLUMNS. raw_delta, log10_raw_delta and
    periodic_delta are on an undocumented convention and are not used downstream.
    """
    rows = _rows(path, "Fig 2")
    header = [str(c).strip() if c is not None else None for c in rows[0]]
    if header[: len(TABLE_COLUMNS)] != list(TABLE_COLUMNS):
        raise ValueError(f"Fig 2: unexpected header {header}")
    body = [r[: len(TABLE_COLUMNS)] for r in rows[1:] if r and r[0] is not None]
    df = pd.DataFrame(body, columns=list(TABLE_COLUMNS.values()))

    if len(df) != N_REC:
        raise ValueError(f"Fig 2: expected {N_REC} rows, got {len(df)}")
    counts = df["genotype"].value_counts().to_dict()
    if counts != {"NT": N_TD, "AS": N_AS}:
        raise ValueError(f"Fig 2: expected {N_TD} NT + {N_AS} AS, got {counts}")
    n_as_ids = df.loc[df["genotype"] == "AS", "id"].nunique()
    if n_as_ids != N_AS_SUBJECTS:
        raise ValueError(f"Fig 2: expected {N_AS_SUBJECTS} unique AS IDs, got {n_as_ids}")
    return df


def map_table_to_spectra(table, reported, tol=1e-5):
    """Map each Fig 2 row to its Fig 1 label by nearest (exponent, offset).

    Never positional: Fig 2 is AS-first, Fig 1 is TD-first, and the order within
    blocks differs. Distance is max(|d exponent|, |d offset|); every row must have
    a match within `tol`, the 344 matches must be unique, and the group (NT->TD,
    AS->AS) must agree for every row. Returns a Series indexed like `table`.
    """
    t = table[["exponent", "offset"]].to_numpy(float)
    r = reported[["exponent", "offset"]].to_numpy(float)
    dist = np.abs(t[:, None, :] - r[None, :, :]).max(axis=2)
    best = dist.argmin(axis=1)
    best_dist = dist[np.arange(len(t)), best]

    far = np.flatnonzero(best_dist > tol)
    if far.size:
        raise ValueError(f"{far.size} Fig 2 rows have no Fig 1 match within {tol} (worst {best_dist.max():.3g})")
    labels = reported.index.to_numpy()[best]
    if len(set(labels)) != len(table):
        raise ValueError(f"Fig 2 -> Fig 1 mapping is not one-to-one ({len(set(labels))} unique labels)")
    group = table["genotype"].map(GENOTYPE_TO_GROUP).to_numpy()
    label_group = np.array([lab[:2] for lab in labels])
    mismatch = np.flatnonzero(group != label_group)
    if mismatch.size:
        raise ValueError(f"{mismatch.size} Fig 2 rows matched a Fig 1 label of the other group")
    return pd.Series(labels, index=table.index, name="label")


def load_all(path):
    """Return (tidy, psd): one row per recording in Fig 1 order, psd[i] <-> tidy.iloc[i].

    tidy columns: label, group ('TD'/'AS'), id, age, sex, del_status, then the
    reported features (REPORTED_FEATURES).
    """
    freqs, psd, labels = load_spectra(path)
    reported = load_reported(path)
    if list(reported.index) != labels:
        raise ValueError("Fig 1a and Fig 1b,d,f,g labels differ")
    table = load_table(path)
    table = table.assign(label=map_table_to_spectra(table, reported))

    demo = table.set_index("label")[["id", "age", "sex", "del_status"]].reindex(labels)
    tidy = demo.join(reported).reset_index()
    tidy.insert(1, "group", tidy["label"].str[:2])
    return tidy, psd


if __name__ == "__main__":
    import sys

    logging.basicConfig(level=logging.INFO)
    tidy, psd = load_all(sys.argv[1] if len(sys.argv) > 1 else "refs/supplementary-data-1.xlsx")
    print(tidy)
    print("psd", psd.shape)
