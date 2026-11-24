# PROTOCOL — predeclared analysis plan (v2.4, 2026-11-19)

Predeclared content froze at v2.3 (2026-11-05). v2.4 adds post-analysis records only, each marked as written 2026-11-19, after all results were seen; no predeclared text was changed.

### Primary estimand

The paired within-recording difference in the integrated 2–4 Hz periodic quantity produced by knee-mode versus fixed-mode specparam fits on the same PSD (fixed minus knee), using the same 1–49.5 Hz fit range. The periodic quantity is the paper's: the band rule (see Conventions) applied to specparam's fitted peak model in log10 units (`config.yaml`: `periodic_quantity_primary`). We will report the median paired difference and a recording-level bootstrap 95% confidence interval of the median, plus paired plots. The spectra are the 344 per-recording PSDs (185 TD, 159 AS) in Supplementary Data 1 of the paper.

### Secondary outcomes

Predeclared: (a) agreement between the knee refit and the paper's reported exponent, offset and periodic delta (validation), measured by Pearson r, mean bias, and 95% limits of agreement, separately for TD and AS; (b) AS-vs-TD Cohen's d for the primary periodic quantity and for exponent under each mode; (c) the TD age association of each feature under each mode (linear first; a GAM only to compare shape with the paper's Fig. 3 TD curve); (d) low-frequency peak counts, fit-error and convergence diagnostics, knee-flag rate (knee ≤ 0 or knee frequency outside 0.5–20 Hz) by group and age bin; and the count of recordings whose peak-model periodic delta is 0 while the raw-residual periodic delta exceeds 0.3, by mode and group (edge-rule drops); (e) the primary paired comparison repeated with the secondary periodic quantity, the band rule applied to the raw log residual log10(PSD) − log10(aperiodic fit) (`periodic_quantity_secondary`). Knee-mode and fixed-mode exponents will not be treated as interchangeable estimates of one parameter.

**(a), continued.** Diagnostic added 2026-11-05, before running: a single knee curve is fitted to log10(Fig 1e) per recording; the maximum residual and the agreement of the fitted (b, χ) with the reported values classify Fig 1e as a single fit on the averaged spectrum or not. The three readings are recorded in NOTES before the run.

### Data and unit of analysis

Supplementary Data 1 sheets used; all 344 recordings; the recording is the unit (as in the paper). The 159 AS recordings come from 95 participants (Fig 2 `ID`); repeated visits are not modeled in the primary analysis, and are handled by the cluster-bootstrap sensitivity analysis under Statistics.

### Inclusion

Primary analysis: all recordings for which both fits converge. Recordings whose knee fit returns knee ≤ 0 or a knee frequency outside 0.5–20 Hz are flagged, counted by group, and retained; the primary comparison is repeated without them as a predeclared sensitivity. QC gate item 2 (“both fits succeed”) refers to convergence, not to the knee flag. No R² threshold.

### Conventions (v2, settled from the data 2026-11-01)

Derived by `src/validate.py` over all 344 recordings; numbers in `config.yaml` (`conventions`), aggregate summary in `tables/delta_rules.csv`. Per-recording derived values are written to `data/` and stay local (CC BY-NC-ND 4.0).

- **Band integration.** Band value = (hi − lo) × mean of the values at the 0.5 Hz bins with lo ≤ f ≤ hi (both endpoints included). This reproduces the paper's reported delta for 344/344 recordings to rounding (max relative error 4.7e-7) and the reported Total 1–50 Hz (= 0.5 × sum of all 98 bins, 1.0–49.5 Hz) for 343/344. The trapezoid rule, the expected candidate, was rejected: median −0.98%, 5th–95th percentile −4.3% to +9.9%, max |error| 18.5%. Theta, alpha, beta and gamma do not follow this rule on closed intervals; they are not used.
- **Delta bins.** 2.0, 2.5, 3.0, 3.5, 4.0 Hz; delta = 2 × mean = 0.4 × sum of these five bins.
- **Paper's periodic delta.** The same rule applied to Fig 1c (the fitted peak model, log10 units) reproduces the reported Fig 1d periodic delta for 344/344 (max |difference| 1.04e-5). It is the fitted peak model, not the raw log residual. This is the primary periodic quantity.
- **Raw log residual (secondary periodic quantity).** log10(PSD) − log10(aperiodic fit), integrated over the delta bins with the same rule. On the paper's own fits (log10 Fig 1a − log10 Fig 1e) it differs from the reported periodic delta (median |difference| 0.091, max 2.23, r = 0.939).
- **Fig 1a vs 1c vs 1e.** log10(Fig 1a) − log10(Fig 1e) − Fig 1c: median |difference| 0.033; per-recording max median 0.125, at a median of 9.5 Hz (IQR 5.5–14.6 Hz). Max in 6–12 Hz for 144/344, 1–4 Hz for 81, 30–50 Hz for 73. Overall max 1.60 at AS154, 2.5 Hz, a sharp delta peak (PSD 1891 µV²/Hz) that the peak model under-fits.
- **Fig 2 → Fig 1 mapping.** Nearest match on (exponent, offset), max-abs distance ≤ 1e-5; never positional. All 344 matched (worst 1.9e-6, TD133); runner-up ≥ 4.4e-3.
- **Data quirks.** TD133 reported offset is `-0.1314*` (trailing asterisk stripped). Fig 2 `Raw_delta` and `Periodic_delta` are on an undocumented convention and are not used.

**Deviations from the paper's reported values:** AS159 Total 1–50 Hz (reported 417.37, ours 416.19, −0.28%); its delta matches exactly. All recorded deviations, including the knee-refit biases from secondary (a), are listed in `tables/deviations.csv` (written by `src/validation.py`).

*Added 2026-11-19.* The same file also lists every difference from this protocol, with its reason; see Deviations from this protocol.

### QC gate

All six must pass before the primary analysis is run.

1. Loader assertions pass.
2. Both fits converge for at least 95% of recordings, with no systematic convergence failure by age or group. Knee flags are counted under item 4 and are not failures.
3. R² distributions in both modes are reported; no threshold is borrowed from the paper.
4. Knee ≤ 0 or knee frequency outside 0.5–20 Hz (knee is unbounded in specparam 2.0.0rc2, so “at bound” cannot occur) is counted by group and age bin and listed.
5. Ten overlays per group (raw PSD, knee and fixed aperiodic fits, both fitted peak models, both raw residuals) are inspected and saved.
6. The whole thing reruns from a clean environment.

### Statistics

Median paired difference; recording-level bootstrap CI of the median, 10,000 resamples, seed 20261029; Welch Cohen's d with bootstrap CI; linear model of each feature on log age per mode, slopes with CIs, slope difference by paired bootstrap.

Sensitivity (predeclared): the primary CI recomputed with a cluster bootstrap that resamples participants (Fig 2 `ID`; 95 AS participants, 185 TD) with all their recordings, same resample count and seed. Reported alongside the recording-level CI; the primary result is the recording-level one.

Sensitivity (predeclared): the primary comparison repeated excluding knee-flagged recordings (knee ≤ 0 or knee frequency outside 0.5–20 Hz).

### Expectations from synthetic validation (2026-11-03, before any real-data analysis)

On the real spectra, fixed-mode periodic delta is expected to exceed knee-mode by roughly 0.3 log10 units in TD, and by more in AS, because AS knee frequencies are lower (paper: 3.2 vs. 4.2 Hz) and exponents steeper. The AS–TD contrast in periodic delta under fixed mode is therefore expected not to vanish and may increase, for a reason unrelated to oscillations.

If the real data do that, the result is a confirmed mechanism, not a post-hoc story. If they don't, that is a finding too, and it will be reported as such.

*Basis.* `tables/sim_summary.csv`, paired fixed − knee periodic delta, medians at noise 0.02 / 0.05: no peak +0.32 / +0.31; exponent steepened by the AS − TD Δχ (knee frequency 3.56 → 2.94 Hz) +0.44 / +0.41; AS-sized 3 Hz peak +0.57 / +0.46; offset, exponent and peak together +0.67 / +0.58. The generator is the knee model with white log10 noise, so these magnitudes are a guide to direction and rough size, not a forecast to the decimal. Condition (c′), which steepens the exponent while holding the knee frequency, produces most of the same fixed-mode excess as (c) (+0.41 vs. +0.44), so the excess comes mainly from the slope rather than the knee position; the expectation above names both because the real AS recordings differ from TD in both.

*2026-11-05, before `analysis.py`.* Knee flags concentrate under age 3 in TD (13 of 17; 11 under 1.1 years). A flagged knee fit coincides with the fixed line, so Δexponent ≈ 0 by construction for those recordings and larger for older children; this may itself produce an age trend in Δexponent. Secondary (c) reports slopes with and without flagged recordings.

*Real-data exposure before this was written.* Knee-mode fits of the 185 TD spectra, used only for their medians of offset, knee and exponent (the synthetic baseline); the TD1 smoke test in both modes. No AS spectrum has been fit, no fixed-mode fit beyond TD1, and no periodic-delta value or AS–TD contrast has been computed from any refit.

2026-11-05. Since the above: knee-mode aggregates per group for exponent, offset and periodic delta (secondary a) and 20 overlays showing both modes per recording. Still unseen: any aggregate of fixed-mode periodic delta, and any fixed − knee difference across recordings or groups.

*Added 2026-11-19.* Nothing predeclared remains unseen: every predeclared outcome, including fixed-mode periodic delta and the fixed − knee differences, was computed on 2026-11-14 (NOTES, that date).

*Outcome of the expectations, written 2026-11-19 after all results were seen* (results of 2026-11-14 and 2026-11-16):
- Fixed exceeds knee by roughly 0.3 in TD: **direction held, size not held.** TD median +0.110 [0.094, 0.129]; the CI excludes 0.3 (`tables/primary.csv`).
- … and by more in AS: **held.** AS median +0.426 [0.349, 0.481].
- The AS–TD contrast under fixed mode does not vanish and may increase: **held; it increased.** Periodic-delta Cohen's d 1.27 [1.08, 1.49] under knee, 1.72 [1.50, 1.98] under fixed (`tables/group_contrast_by_mode.csv`).
- Because AS knee frequencies are lower and exponents steeper: **not tested, so the mechanism is not confirmed.** The results measure the change in the contrast, not its cause.
- 2026-11-05, Δexponent ≈ 0 for knee-flagged recordings: **held.** TD +0.028 [−0.008, +0.067] (n 17), AS −0.000 [−0.039, +0.050] (n 8) (`tables/secondary_b.csv`).
- 2026-11-05, the flagged recordings may themselves produce an age trend in Δexponent: the trend is present (TD −0.241 per doubling of age [−0.274, −0.209]) but **not produced by them**: −0.227 [−0.266, −0.191] without them (`tables/secondary_c.csv`).
- Why the TD size was overpredicted, **post hoc** (NOTES 2026-11-16; `tables/sim_gap_twins.csv`): the synthetic generator lacked each recording's own peaks (45% of the 0.21 difference) and the structure of the real residual (32%); using one median curve (15%) and the residual's size (8%) mattered less.

### Exploratory, not predeclared

Fixed mode 2–40 Hz; doublexp; the Fig 6 mouse PCA; HBN.

Exploratory analyses, if reported, are labeled as such and fall outside the claims below.

*Added 2026-11-19.* None of the four analyses listed above was run. These were run instead; all are labelled in their tables and fall outside the claims below:
- exploratory (2026-11-14): AS age slopes under each mode; the AS − TD median periodic delta under each mode and its fixed − knee change; the fixed − knee change in exponent and offset by group, with knee-flagged and unflagged subsets; the mean paired difference (descriptive);
- post hoc (2026-11-14): the decomposition of the gap between refit and reported periodic delta, its subsets, and cross-tabs of the Fig 1e residual pattern against the knee flag and age under 3;
- post hoc (2026-11-16): why the synthetic validation overpredicted the TD size, using synthetic twins of each TD recording.

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

### Deviations from this protocol

*Added 2026-11-19, after the analysis.* Every difference between what was done and this protocol's predeclared content (as of v2.3) is listed in `tables/deviations.csv`, in the rows whose `kind` is not `reported value`, with the reason recorded at the time or a note that none was. Those that touch predeclared items:
- Secondary (c): the GAM shape comparison with the paper's Fig. 3 TD curve was not run, and no reason was recorded. Whether the paper's TD curve keeps its shape under fixed mode is answered for the linear slope only.
- Primary, paired plots: the committed figure bins recordings by pair mean and shows the median and IQR of the paired difference per bin; the per-recording paired plot stays local, because per-recording derived values are not redistributed (CC BY-NC-ND 4.0).
- Secondary (b) is reported in `tables/group_contrast_by_mode.csv`; `tables/secondary_b.csv` holds the exploratory fixed − knee parameter changes.
- Statistics, "log age": log2(age), so slopes are per doubling of age; the base was not fixed here.
- QC gate item 2 failed as written on 2026-11-04 and was resolved by the v2.3 rewording to convergence (344/344 in both modes).

The estimand, the inclusion rule, the bootstrap and both sensitivity analyses were carried out as written.

### Change log

- v1 2026-10-31: initial.
- v2 2026-11-01: conventions filled from the data. Expected the trapezoid rule (off by 1–5%); the data show the exact rule (hi − lo) × mean over closed-interval bins, so trapezoid was rejected. Found that the paper's periodic delta is the band rule over the fitted peak model (Fig 1c), not the raw log residual. Estimand clarified: the primary periodic quantity is the paper's (band rule over the fitted peak model), "summed" became "integrated" (the band rule is 2 × mean, not a sum), and the bootstrap is recording-level, not participant-level; the raw log residual moved to secondary (e); participant cluster-bootstrap sensitivity added. QC overlays now include both fitted peak models; secondary (a) states its agreement measures; secondary (b) names the primary periodic quantity. Recorded AS159 Total as a deviation. Loader with assertions (`src/load_supplement.py`, `tests/test_loader.py`).
- v2.1 2026-11-03: added "Expectations from synthetic validation", written after the synthetic grid and before any real-data analysis. Estimand, outcomes, QC gate and statistics unchanged.
- 2026-11-04: pointer from the deviations sentence to `tables/deviations.csv`. Bookkeeping only; nothing else changed.
- v2.2 2026-11-04: inclusion rule amended before any AS fit. The TD dry run flagged 17/185 knee fits (9 knee ≤ 0, 8 knee frequency outside 0.5–20 Hz); these are not convergence failures, and excluding them conditions on one arm of the comparison. Flagged recordings are now retained in the primary analysis and dropped only in a predeclared sensitivity. QC gate item 2 clarified to mean convergence.
- v2.3 2026-11-05: QC item 2 reworded to convergence (it was the v2.2 intent; the text still said “succeed”); item 4 and secondary (d) reworded to the knee flag, since the knee is unbounded; edge-rule count added to (d); Fig 1e single-curve diagnostic added to (a); expectation on the age pattern of knee flags added. Estimand and statistics unchanged. v2.2 was decided 2026-11-04 before the AS fits were run and is being committed 2026-11-05 after them; the decision predates the data.
- v2.4 2026-11-19: post-analysis records, written after all results were seen (analysis 2026-11-14, follow-ups 2026-11-16). Added, each marked as written 2026-11-19: a note under the title; a pointer from the deviations sentence to the protocol rows of `tables/deviations.csv`; a closing line to the real-data exposure note; the outcome of each expectation; the exploratory and post hoc analyses that were run; and the section Deviations from this protocol. No predeclared text was changed: the estimand, secondary outcomes (including the GAM clause in (c), which was not run and is recorded as a deviation, not removed), inclusion, conventions, QC gate, statistics, expectations and the claims section read as in v2.3.
