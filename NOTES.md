# NOTES

## 2026-10-26 — Supplement: hash, license, workbook traps

*Reconstructed 2026-11-01 from the repo and the loader spec; details observed on the day not recorded.*

- File: Supplementary Data 1 (`refs/supplementary-data-1.xlsx`), SHA-256 `aa6dfb4cdd2a47b98cf2401f15a4ac609104eb8b45df6508c3b54519c37b8736` (re-checked 2026-11-01). The local copy's modification time is 2026-10-29 17:53, so it was re-saved or re-downloaded after this entry.
- License: CC BY-NC-ND 4.0 (from the paper). Consequence: per-recording derived values stay local (`data/`, gitignored); only aggregates are published (`tables/`).
- Traps found Oct 25–26, each now a loader assertion: unlabeled columns after the AS block in Fig 1a; a trailing `*` on one reported value (TD133 offset); Fig 2 is AS-first while Fig 1 is TD-first, with a different order within blocks, so rows cannot be matched by position.

## 2026-10-29 — Smoke test: refit TD1 with the paper's settings

*Reconstructed 2026-11-01; numbers are from a rerun of `smoke_test3.py` on 2026-11-01 (specparam 2.0.0rc2), not the Oct 29 output.*

`smoke_test.py` (fooof 1.1) → `smoke_test2.py` (specparam port) → `smoke_test3.py` (pinned to specparam 2.0.0rc2, explicit targets, knee frequency = k^(1/exponent)). TD1, fit range 1–49.5 Hz, peak_width_limits 0.5–12, peak_threshold 2.0, max_n_peaks inf, min_peak_height 0:

| | reported | knee refit | fixed refit |
|---|---|---|---|
| exponent | 1.9809 | 1.9350 | 1.5288 |
| offset | 2.3131 | 2.3077 | 1.7351 |
| R² | | 0.9985 | 0.9906 |
| peaks | | 13 | 3 |

Knee frequency 2.84 Hz (k = 7.54). All four targets pass, but the targets bracket the refit, not the paper: the knee refit's exponent is 0.046 below the reported value. That refit-vs-reported gap is what secondary (a) measures across all 344.

## 2026-11-01 — Loader with assertions; band conventions settled

`src/load_supplement.py`, `tests/test_loader.py` (13 tests, all pass).

- Fig 1a: 98 frequencies 1.0–49.5 Hz step 0.5; labels TD1…TD185 then AS1…AS159; 26 empty trailing columns skipped; all PSD values finite and > 0 (min 0.0077).
- Fig 1b,d,f,g: one starred cell, TD133 offset `-0.1314*`; stripped and logged.
- Fig 2: 344 rows, 185 NT + 159 AS, 95 unique AS IDs.
- Fig 2 → Fig 1 mapping, nearest on (exponent, offset), max-abs distance: 344/344 matched; worst 1.9e-6 (TD133, whose offset is rounded to 4 decimals); runner-up for every row ≥ 4.4e-3. The mapping is not a block swap: 122 TD rows are reordered (from Fig 2 row 222 on).

Band rule (`src/validate.py`; `tables/delta_rules.csv`). Delta vs reported Fig 1b, relative error = ours / reported − 1:

| rule | median | max \|e\| | closest for |
|---|---|---|---|
| 2·mean[2,4] = 0.4·Σ | 2e-16 | 4.7e-7 | 344 |
| trapezoid[2,4] | −0.98% | 18.5% | 0 |
| 0.5·Σ[2,4) | +7.8% | 22% | 0 |
| 0.5·Σ[2,4] | +25.0% (every recording) | 25% | 0 |

Expected the trapezoid, off by 1–5%; the data show an exact rule, (hi − lo) × mean over the closed-interval bins. Total 1–50 Hz = 0.5·Σ (same rule, 49 Hz × mean) for 343/344; AS159 −0.28% (reported 417.37, ours 416.19), recorded as a deviation.

## 2026-11-01 — The paper's "periodic delta" is the fitted peak model, not the residual

**What the paper says.** Methods, "EEG power spectral analysis":

> "For each power spectrum, we extracted the aperiodic component using specparam and subtracted this component from the raw PSD to obtain the periodic component. Because the subtraction was performed in log space, the resulting periodic component reflects a power ratio rather than an absolute value, and is therefore unitless. Periodic power within each frequency band was calculated by summing within the same bands as above."

The mouse section says the same: "The periodic power spectrum was calculated as log10(PSD)—log10(aperiodic power)." The project's v1 periodic quantity followed this text: the raw log residual, log10(PSD) − log10(aperiodic fit), integrated over 2–4 Hz.

**The check.** With the band rule already pinned down by the reported delta (2 × mean of the 0.5 Hz bins 2.0, 2.5, 3.0, 3.5, 4.0 Hz; 344/344 exact to 4.7e-7), apply the same rule to two candidates and compare each with the reported Fig 1d periodic delta, all 344 recordings (`src/validate.py`, check 2b):

```
periodic_delta_hat = 2 * mean(x[:, 2.0 <= f <= 4.0])

x = Fig 1c                          median |d| 4.0e-07   max |d| 1.04e-05   r = 1.0000
x = log10(Fig 1a) - log10(Fig 1e)   median |d| 9.1e-02   max |d| 2.23       r = 0.9389
```

The supplement's Fig 1c is not that residual; it is the fitted peak model, and the paper's band rule over it reproduces Fig 1d to 1e-5. The raw residual on the paper's own fits does not (r = 0.94, median |diff| 0.09). The project follows the source data and records this as a documentation note.

The two candidates differ by the model misfit: |log10(1a) − log10(1e) − 1c| has median 0.033, and its largest value, 1.60, is at AS154, 2.5 Hz, a sharp delta peak (PSD 1891 µV²/Hz) that the Gaussian peak model under-fits.

**Two text-vs-data differences, different in kind.**

- *Raw residual vs. fitted peak model: a difference in kind.* Only the fitted Gaussian peaks reproduce the reported periodic delta. specparam's own documentation also calls the peak fit "the periodic component", so the authors' usage is defensible shorthand; the Methods sentence is a simplification, not a mistake.
- *"Summing" vs. (hi − lo) × mean: a difference of a constant.* For delta the rule is 0.4 × Σ, a scaled sum, so "summing" is true up to the factor; the factor still has to be matched to reproduce the numbers.

**Why it matters.** The estimand changes. Under the paper's quantity, "periodic delta" is whatever the peak fitter decides is a peak; a knee-vs-fixed difference then runs through two fitted components (aperiodic and peaks), not one. The raw residual would have measured something the paper never reported. PROTOCOL v2 makes the paper's quantity primary and the raw residual secondary (e).

**Essay.** Best "one honest difficulty" candidate so far: the text and the data differed, and it was resolved by checking. The difference is invisible in the summary statistics (r = 0.94) and only showed up because the supplement included the fitted peak model as its own sheet.

**Email to the authors (decide at the freeze).** One question about how the reported features were computed from the region-averaged PSDs, in two parts: the periodic definition, and the AS exponent gap.
