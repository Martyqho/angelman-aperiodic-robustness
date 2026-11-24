# angelman-aperiodic-robustness

How much does the choice of aperiodic model (knee vs fixed) change the
"periodic delta" EEG measure reported for Angelman syndrome?

**Status:** Phases 1–4 complete (loader, synthetic validation, full refit
of 344 spectra, QC gate 6/6, primary analysis). Analysis frozen at tag v1.0-freeze (2026-11-16 results; tagged 2026-11-24);
results below.
Analysis plan: [PROTOCOL.md](PROTOCOL.md). Lab notebook: [NOTES.md](NOTES.md).

## What this is

A methodological robustness study, not a clinical replication. The target
paper is

> Montgomery, Shide, Didouchevski, Dickinson & Sidorov (2026). *Periodic
> and aperiodic contributions to EEG delta power are translatable and
> complementary Angelman syndrome biomarkers.* Communications Medicine
> 6:532. https://www.nature.com/articles/s43856-026-01741-0
> (preprint: https://doi.org/10.1101/2025.09.08.674941)

The paper separates each recording's power spectrum into periodic and
aperiodic parts with specparam (v2.0.0rc2) in **knee** mode and reports
"periodic delta" (2–4 Hz) as a biomarker. This project refits the paper's
own 344 per-recording spectra (185 TD, 159 AS, from Supplementary Data 1)
in both knee and fixed mode and asks one predeclared question:

> **Primary estimand.** The within-recording difference (fixed − knee) in
> the paper's periodic-delta measure; median across recordings with a
> recording-level bootstrap 95% CI (10,000 resamples, seed 20261029),
> overall and by group.

Secondary analyses, sensitivities (knee-flagged recordings excluded;
cluster bootstrap by participant where possible), expectations written
before the data were seen, and a change log are in PROTOCOL.md.

## Data

The spectra are **Supplementary Data 1** of the paper, downloaded from the
"Supplementary information" section of the article page above. The file is
licensed CC BY-NC-ND 4.0 and is **not redistributed here**. Code is MIT-licensed; see License below. To run the
pipeline, download it to `refs/supplementary-data-1.xlsx`. The pipeline
verifies it before anything else:

    SHA-256  aa6dfb4cdd2a47b98cf2401f15a4ac609104eb8b45df6508c3b54519c37b8736
    size     2,508,461 bytes

Per-recording derived values (fits, peaks, features) are written to `data/`
and stay local (gitignored). Only aggregate tables and figures are
committed, in `tables/` and `figures/`.

Conventions reverse-engineered from the supplement and documented in
`config.yaml` and NOTES.md: band power = (hi − lo) × mean over the closed
0.5 Hz bins (exact for 344/344 delta values); the paper's "periodic"
spectrum (Fig 1c) is specparam's fitted peak model, not the raw residual.

## Running it

    python -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt        # pinned; specparam 2.0.0rc2 is the version the paper used
    bash scripts/run_all.sh                # sha check → validate → fit_spectra → qc → analysis → validation
                                           # → fig1e_diagnostic → gap_subsets → sim_gap_twins → pytest
    bash scripts/run_all.sh --with-sim     # also reruns the synthetic validation
Tested with Python 3.13 on macOS; `requirements.txt` is a full `pip freeze` from that environment.

Settings (fit range 1–49.5 Hz, peak_width_limits 0.5–12, peak_threshold
2.0, max_n_peaks inf, min_peak_height 0) are the paper's and live in
`config.yaml`; nothing is hard-coded in the scripts.

## Layout

    LICENSE              MIT, for the code
    PROTOCOL.md          predeclared analysis plan, versioned, with change log
    NOTES.md             dated lab notebook: decisions, findings, retractions
    config.yaml          paper citation, supplement hash, grid, bands, specparam settings, conventions
    scripts/run_all.sh   end-to-end pipeline
    src/
      config.py            loads config.yaml
      load_supplement.py   reads the xlsx, maps Fig 2 rows to Fig 1 by value, never by position
      validate.py          band-rule and convention checks against reported values
      fit_spectra.py       knee and fixed fits for every recording (fit_one)
      qc.py                convergence, knee flags, overlays, fit-quality tables
      analysis.py          primary, sensitivities, secondaries (a)–(e), group contrast by mode
      validation.py        refit vs reported exponent, offset, periodic delta; tables/deviations.csv
      simulate.py          synthetic spectra with known parameters (Phase 2)
      plot_simulation.py   figures/fig_synthetic_validation.png
      fig1e_diagnostic.py  is the paper's aperiodic component a single knee curve?
      gap_subsets.py       periodic-gap subsets and Fig 1e cross-tabs (post hoc)
      sim_gap_twins.py     why the synthetic validation overpredicted the TD gap (post hoc)
    tests/               pytest; run by run_all.sh
    tables/              aggregate results (committed)
    figures/             committed figures
    data/, refs/         local only

## Results

Every number below is from a table in `tables/`;
the Phase 4 tables label each row predeclared, exploratory or post hoc. Periodic delta
is the paper's quantity (band rule over specparam's fitted peak model, 2–4 Hz), in log10
units.

### Primary

d = periodic delta under fixed mode minus under knee mode, per recording; median with a
recording-level bootstrap 95% CI (10,000 resamples, seed 20261029). `tables/primary.csv`.

| | n | median d [95% CI] | recordings with d > 0 |
|---|---|---|---|
| All | 344 | **+0.197 [0.161, 0.229]** | 92% |
| TD | 185 | +0.110 [0.094, 0.129] | 89% |
| AS | 159 | +0.426 [0.349, 0.481] | 96% |

Fitting the same spectra in fixed instead of knee mode raises periodic delta in almost
every recording, about four times as much in AS as in TD.

### Sensitivity (predeclared)

- **Without the 25 knee-flagged recordings** (17 TD, 8 AS): all +0.217 [0.183, 0.241],
  TD +0.123 [0.102, 0.144], AS +0.446 [0.394, 0.502].
- **Participant cluster bootstrap** (95 AS participants; TD has one recording each, so its
  CI is unchanged by construction): all [0.158, 0.233], AS [0.331, 0.482].

Neither changes the conclusion.

### AS–TD contrast under each mode

`tables/group_contrast_by_mode.csv`.

| | knee mode | fixed mode |
|---|---|---|
| Periodic delta, Welch Cohen's d (predeclared) | 1.27 [1.08, 1.49] | 1.72 [1.50, 1.98] |
| Exponent, Welch Cohen's d (predeclared) | 0.87 [0.66, 1.11] | 1.53 [1.32, 1.76] |
| Periodic delta, median AS − TD (exploratory) | +0.47 [0.38, 0.65] | +0.72 [0.62, 0.89] |

The contrast is not stable under the change of model: it is larger in fixed mode for both
features (median difference up by +0.24 [0.14, 0.34], exploratory). This measures how much
the contrast depends on the model; it does not say which model is right.

### Fig 1e

None of the 344 aperiodic curves the paper publishes (Fig 1e) is a single specparam knee
curve: the closest single knee curve misses every one by far more than spreadsheet rounding
(at least 2.0e-4 log10 units, median 7.4e-3, largest at the lowest frequencies). The reported features
therefore come from some computation other than one knee fit of the published
region-averaged spectrum, which this study cannot identify; our knee refit matches the
reported exponent on average (bias −0.012 TD, −0.002 AS) with per-recording limits of
agreement of about ±0.2 (TD) and ±0.35 (AS), and gives AS periodic delta 0.125 lower, a gap
that is reported here, not explained.

### Further detail

Secondary outcomes (age slopes, diagnostics), the scorecard of the expectations written
before the data were seen, and the follow-up on why the synthetic validation overpredicted
the TD size are in [NOTES.md](NOTES.md) (2026-11-14, 2026-11-16). Every difference from
PROTOCOL, with its reason where one was recorded, is in `tables/deviations.csv`.

### What this study cannot claim

Copied from [PROTOCOL.md](PROTOCOL.md) (v2.3; unchanged in v2.4).

This is a methodological robustness study on the authors' published spectra. Its scope is the step from a
per-recording PSD to the reported spectral features, and nothing else.

It **cannot** say:
- anything about the clinical validity of delta power or aperiodic features as Angelman syndrome biomarkers;
- anything about processing upstream of the PSD: filtering, artifact rejection, channel selection, region
  averaging, or the order in which channels were averaged and fit. Those are the authors' and are taken as given;
  the refit-vs-reported gap is reported, not explained;
- anything about the Bayley associations, the genotype analysis, the mouse recordings, or behavior;
- which fit mode is "correct". Knee- and fixed-mode exponents are different quantities by construction, so a
  difference between them is a statement about model dependence, not an error in either;
- the true developmental trajectory of the aperiodic exponent in typically developing children; only whether
  the paper's TD curve keeps its shape under the other model on this cohort;
- anything about recordings as individuals. The recording is the unit, repeated AS visits are not modeled as
  such in the primary analysis, and results are reported as aggregates;
- generalization beyond this cohort, these channels, this frequency range (1–49.5 Hz) and this specparam
  version (2.0.0rc2).

## License

The code is released under the MIT License; see [LICENSE](LICENSE). That covers `src/`,
`tests/`, `scripts/`, `config.yaml`, `pyproject.toml` and `requirements.txt`, and nothing
else. Supplementary Data 1 belongs to the paper's authors, is licensed CC BY-NC-ND 4.0 and
is not included here. The documents (README, PROTOCOL, NOTES), tables and figures are © 2026 Marty Ho. They may be quoted and cited with attribution; they are not under the MIT License because the tables and figures are derived from the CC BY-NC-ND supplement.

## Author

Marty Ho, Alpharetta/GA, Innovation Academy junior. Independent project,
2026. Questions: marty.q.ho@gmail.com
