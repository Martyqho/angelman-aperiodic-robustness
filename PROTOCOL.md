# PROTOCOL — predeclared analysis plan (v2, 2026-11-01)

### Primary estimand

The paired within-recording difference in the integrated 2–4 Hz periodic quantity produced by knee-mode versus fixed-mode specparam fits on the same PSD (fixed minus knee), using the same 1–49.5 Hz fit range. The periodic quantity is the paper's: the band rule (see Conventions) applied to specparam's fitted peak model in log10 units (`config.yaml`: `periodic_quantity_primary`). We will report the median paired difference and a recording-level bootstrap 95% confidence interval of the median, plus paired plots. The spectra are the 344 per-recording PSDs (185 TD, 159 AS) in Supplementary Data 1 of the paper.

### Secondary outcomes

Predeclared: (a) agreement between the knee refit and the paper's reported exponent, offset and periodic delta (validation), measured by Pearson r, mean bias, and 95% limits of agreement, separately for TD and AS; (b) AS-vs-TD Cohen's d for the primary periodic quantity and for exponent under each mode; (c) the TD age association of each feature under each mode (linear first; a GAM only to compare shape with the paper's Fig. 3 TD curve); (d) low-frequency peak counts, fit-error and convergence diagnostics, knee-at-bound rate; (e) the primary paired comparison repeated with the secondary periodic quantity, the band rule applied to the raw log residual log10(PSD) − log10(aperiodic fit) (`periodic_quantity_secondary`). Knee-mode and fixed-mode exponents will not be treated as interchangeable estimates of one parameter.

### Data and unit of analysis

Supplementary Data 1 sheets used; all 344 recordings; the recording is the unit (as in the paper). The 159 AS recordings come from 95 participants (Fig 2 `ID`); repeated visits are not modeled in the primary analysis, and are handled by the cluster-bootstrap sensitivity analysis under Statistics.

### Inclusion

A recording is excluded from the paired comparison only if a fit fails (no convergence, or knee frequency outside 0.5–20 Hz); excluded recordings are counted and listed. No R² threshold.

### Conventions (v2, settled from the data 2026-11-01)

Derived by `src/validate.py` over all 344 recordings; numbers in `config.yaml` (`conventions`), aggregate summary in `tables/delta_rules.csv`. Per-recording derived values are written to `data/` and stay local (CC BY-NC-ND 4.0).

- **Band integration.** Band value = (hi − lo) × mean of the values at the 0.5 Hz bins with lo ≤ f ≤ hi (both endpoints included). This reproduces the paper's reported delta for 344/344 recordings to rounding (max relative error 4.7e-7) and the reported Total 1–50 Hz (= 0.5 × sum of all 98 bins, 1.0–49.5 Hz) for 343/344. The trapezoid rule, the expected candidate, was rejected: median −0.98%, 5th–95th percentile −4.3% to +9.9%, max |error| 18.5%. Theta, alpha, beta and gamma do not follow this rule on closed intervals; they are not used.
- **Delta bins.** 2.0, 2.5, 3.0, 3.5, 4.0 Hz; delta = 2 × mean = 0.4 × sum of these five bins.
- **Paper's periodic delta.** The same rule applied to Fig 1c (the fitted peak model, log10 units) reproduces the reported Fig 1d periodic delta for 344/344 (max |difference| 1.04e-5). It is the fitted peak model, not the raw log residual. This is the primary periodic quantity.
- **Raw log residual (secondary periodic quantity).** log10(PSD) − log10(aperiodic fit), integrated over the delta bins with the same rule. On the paper's own fits (log10 Fig 1a − log10 Fig 1e) it differs from the reported periodic delta (median |difference| 0.091, max 2.23, r = 0.939).
- **Fig 1a vs 1c vs 1e.** log10(Fig 1a) − log10(Fig 1e) − Fig 1c: median |difference| 0.033; per-recording max median 0.125, at a median of 9.5 Hz (IQR 5.5–14.6 Hz). Max in 6–12 Hz for 144/344, 1–4 Hz for 81, 30–50 Hz for 73. Overall max 1.60 at AS154, 2.5 Hz, a sharp delta peak (PSD 1891 µV²/Hz) that the peak model under-fits.
- **Fig 2 → Fig 1 mapping.** Nearest match on (exponent, offset), max-abs distance ≤ 1e-5; never positional. All 344 matched (worst 1.9e-6, TD133); runner-up ≥ 4.4e-3.
- **Data quirks.** TD133 reported offset is `-0.1314*` (trailing asterisk stripped). Fig 2 `Raw_delta` and `Periodic_delta` are on an undocumented convention and are not used.

**Deviations from the paper's reported values:** AS159 Total 1–50 Hz (reported 417.37, ours 416.19, −0.28%); its delta matches exactly.

### QC gate

All six must pass before the primary analysis is run.

1. Loader assertions pass.
2. Both fits succeed for at least 95% of recordings, with no systematic failure by age or group.
3. R² distributions in both modes are reported; no threshold is borrowed from the paper.
4. Knee-at-bound or implausible knee frequency (< 0.5 Hz or > 20 Hz) is counted and listed.
5. Ten overlays per group (raw PSD, knee and fixed aperiodic fits, both fitted peak models, both raw residuals) are inspected and saved.
6. The whole thing reruns from a clean environment.

### Statistics

Median paired difference; recording-level bootstrap CI of the median, 10,000 resamples, seed 20261029; Welch Cohen's d with bootstrap CI; linear model of each feature on log age per mode, slopes with CIs, slope difference by paired bootstrap.

Sensitivity (predeclared): the primary CI recomputed with a cluster bootstrap that resamples participants (Fig 2 `ID`; 95 AS participants, 185 TD) with all their recordings, same resample count and seed. Reported alongside the recording-level CI; the primary result is the recording-level one.

### Exploratory, not predeclared

Fixed mode 2–40 Hz; doublexp; the Fig 6 mouse PCA; HBN.

Exploratory analyses, if reported, are labeled as such and fall outside the claims below.

### What this study cannot claim

This is a methodological robustness study on the authors' published spectra. Its scope is the step from a
per-recording PSD to the reported spectral features, and nothing else.

It **can** say:
- whether the paper's periodic delta (band rule over specparam's fitted peak model) and its aperiodic
  parameters change when the same 344 PSDs are fit in fixed instead of knee mode, and by how much;
- whether the published AS–TD contrasts and the TD age association for those features are stable under
  that change of model;
- how closely a refit with the paper's settings reproduces the paper's reported values, and where it does not.

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

Where the paper's text and its source data disagree (the periodic component is described as the raw residual
but reported as the fitted peak model; band values follow a closed-interval mean rule), this study follows the
source data and records the difference as a documentation note, not as an error.

### Change log

- v1 2026-10-31: initial.
- v2 2026-11-01: conventions filled from the data. Expected the trapezoid rule (off by 1–5%); the data show the exact rule (hi − lo) × mean over closed-interval bins, so trapezoid was rejected. Found that the paper's periodic delta is the band rule over the fitted peak model (Fig 1c), not the raw log residual. Estimand clarified: the primary periodic quantity is the paper's (band rule over the fitted peak model), "summed" became "integrated" (the band rule is 2 × mean, not a sum), and the bootstrap is recording-level, not participant-level; the raw log residual moved to secondary (e); participant cluster-bootstrap sensitivity added. QC overlays now include both fitted peak models; secondary (a) states its agreement measures; secondary (b) names the primary periodic quantity. Recorded AS159 Total as a deviation. Loader with assertions (`src/load_supplement.py`, `tests/test_loader.py`).
