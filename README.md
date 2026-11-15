# angelman-aperiodic-robustness

How much does the choice of aperiodic model (knee vs fixed) change the
"periodic delta" EEG measure reported for Angelman syndrome?

**Status:** Phases 1–3 complete (loader, synthetic validation, full refit
of 344 spectra, QC gate 6/6). Phase 4 (primary analysis) in progress.
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
licensed CC BY-NC-ND 4.0 and is **not redistributed here**. To run the
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
    bash scripts/run_all.sh                # sha check → validate → fit_spectra → qc → validation → pytest
    bash scripts/run_all.sh --with-sim     # also reruns the synthetic validation
Tested with Python 3.13 on macOS; `requirements.txt` is a full `pip freeze` from that environment.

Settings (fit range 1–49.5 Hz, peak_width_limits 0.5–12, peak_threshold
2.0, max_n_peaks inf, min_peak_height 0) are the paper's and live in
`config.yaml`; nothing is hard-coded in the scripts.

## Layout

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
      validation.py        refit vs reported exponent, offset, periodic delta
      simulate.py          synthetic spectra with known parameters (Phase 2)
      plot_simulation.py   figures/fig_synthetic_validation.png
      fig1e_diagnostic.py  is the paper's aperiodic component a single knee curve?
    tests/               pytest; run by run_all.sh
    tables/              aggregate results (committed)
    figures/             committed figures
    data/, refs/         local only

## Results

*To be filled at the analysis freeze.*

## What this study cannot claim

- It is not a replication of the paper's clinical findings and makes no
  claim about biomarker validity, group separation, or associations with
  developmental scores.
- It cannot say which aperiodic mode is "correct"; it measures how much the
  choice matters for one quantity.
- The supplement contains only region-averaged spectra. The paper's
  reported per-recording values are averages of regional fits (see NOTES,
  Fig 1e diagnostic), so exact agreement with them is neither expected nor
  tested; all conclusions concern fits of the region-averaged spectra.
- Nothing here concerns the mouse data or any individual participant.

## Author

Marty Ho, Alpharetta/GA, Innovation Academy junior. Independent project,
2026. Questions: marty.q.ho@gmail.com
