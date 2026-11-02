# Knee vs. fixed aperiodic fits on the Angelman EEG spectra
Montgomery et al. (2026)

A methodological robustness study of the spectral-parameterization pipeline in the paper authored by 
Montgomery, Shide, Didouchevski, Dickinson & Sidorov, *Communications
Medicine* 6:532 (2026),
using the per-recording spectra the authors published as Supplementary Data 1.

## Question

The paper splits elevated
2–4 Hz delta power in
Angelman syndrome into a periodic (oscillatory) part and an aperiodic
(1/f-like) part by fitting
specparam in "knee" mode. How much do the periodic delta residual and
the group contrasts change
when the same spectra are fit in "fixed" mode over the same frequency range?

## Primary result

Pending for now

Across 344 recordings, the summed 2–4
Hz periodic residual was
smaller under fixed mode

![Primary figure](figures/fig_primary_knee_vs_fixed.png)

## What this is and is not

- It is: a refit of the authors' published spectra under two aperiodic
models, a validation of
  the refit against the authors' reported values, and a synthetic-data
check of the method.
- It is not: a replication of the clinical study, of anything upstream
of the Power spectral densities (PSD) filtering,
  artifact rejection, channel averaging, of the Bayley analysis, or
of the mouse work.
- Predeclared protocol: [PROTOCOL.md](PROTOCOL.md).

## Data

The spectra are **not** in this repository. Download Supplementary
Data 1 from the article page
(https://www.nature.com/articles/s43856-026-01741-0, "Supplementary
information" section) and
save it as `refs/supplementary-data-1.xlsx`. Expected SHA-256:

    aa6dfb4cdd2a47b98cf2401f15a4ac609104eb8b45df6508c3b54519c37b8736
(2,508,461 bytes)

The article and its supplementary files are licensed CC BY-NC-ND 4.0.
This repository publishes
aggregate results and the code that regenerates everything from the
authors' file.

## How to run

    # run all the tests
    python3 -m pytest tests/                              
    # smoke test
    python3 src/smoke_test.py [path-to-data]
→ fits → validation → analysis → figures

Runs in a few minutes on a laptop. Python 3.12; specparam 2.0.0rc2;
versions in `requirements.txt`.

## Layout

| Path | Purpose |
| --- | --- |
| `src/load_supplement.py` | Reads the workbook; asserts the frequency grid, the 344 recordings, and the Fig 2 ↔ Fig 1 mapping |
| `src/simulate.py` | Synthetic spectra with known ground truth |
| `src/fit_spectra.py` | specparam knee and fixed fits on all spectra |
| `src/validate.py` | Refit vs. reported values; band-power conventions |
| `src/analysis.py` | Paired differences, bootstrap CIs, effect sizes by mode, age models |
| `tables/`, `figures/` | Aggregate outputs (committed) |
| `NOTES.md` | Dated lab notebook |

## Citation

Montgomery DP, Shide JJ, Didouchevski AM, Dickinson AH, Sidorov MS.
Periodic and aperiodic
contributions to EEG delta power are translatable and complementary
Angelman syndrome biomarkers.
*Commun Med* 6, 532 (2026). https://doi.org/10.1038/s43856-026-01741-0

Donoghue T, et al. Parameterizing neural power spectra into periodic
and aperiodic components.
*Nat Neurosci* 23, 1655–1665 (2020). Software:
https://github.com/fooof-tools/fooof

## License and acknowledgements

Code: MIT (see LICENSE). Data: the authors', under CC BY-NC-ND 4.0,
not redistributed here.

Written by Marty, a high school student, as an independent project.
Methodological discussion,
code review and debugging help from a parent and AI assistants; all
analysis code and text are
the author's.
