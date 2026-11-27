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

**Email to the authors (decide at the freeze).** One question about how the reported features were computed from the region-averaged PSDs, in two parts: the periodic definition, and the AS exponent gap. *(2026-11-04: the AS exponent gap was checked across all 159 AS recordings and not found as a bias; see the validation entry. The AS periodic-delta discrepancy found there is the stronger candidate for the second part.)*

## 2026-11-03 — Synthetic validation: the peak-model quantity separates peaks from aperiodic shifts (in knee mode)

`src/fit_spectra.py` (`fit_one`, the only fit call), `src/simulate.py`, `src/plot_simulation.py`, `tests/test_simulate.py`; `figures/fig_synthetic_validation.png`, `tables/sim_summary.csv`. Per-spectrum results in `data/sim_spectra.csv` (local; regenerate with `python src/simulate.py`, ~12 min).

**fit_one components.** specparam 2.0.0rc2 exposes `fm.get_model('aperiodic'|'peak'|'full', 'log')` (backed by `_ap_fit`, `_peak_fit`, `modeled_spectrum_`). Checked on every fit: log10(psd) − ap_fit == `_spectrum_flat`, and ap_fit + peak_fit == full model.

**Baseline (dry run of the Phase 3 loop).** Knee-mode fit_one over the 185 TD spectra: 185/185 succeed; 9 have knee ≤ 0 (no knee frequency); 8 knee frequencies fall outside 0.5–20 Hz. Medians: b₀ = 2.3597, k₀ = 16.519, χ₀ = 2.2088, knee frequency 3.62 Hz (n = 176; k₀^(1/χ₀) = 3.56 Hz). Paper's TD mean is 4.15 Hz.

**Shifts from the reported values.** Δb = +0.842 and Δχ = +0.393 (AS − NT medians of the Fig 2 offset and exponent). Reported periodic delta medians AS 0.663, TD 0.064 → injected band-rule value 0.599 → A = 0.603 for a Gaussian at 3 Hz, σ = 0.5 (band rule of a unit-height peak = 0.9935). A/2 and 2A also run.

**Generator.** log10 P = b − log10(k + f^χ) + A·exp(−(f − 3)²/(2·0.5²)) + ε, ε ~ Normal(0, nlv) in log10, nlv ∈ {0.02, 0.05}; white noise, not the chi-square noise of a real Welch estimate. Identical (max |diff| 0.0) to `specparam.sim.sim_power_spectrum` with the same parameters and seed, with and without noise: log10 and σ-as-standard-deviation conventions agree. Noise-free: knee fit returns b₀, k₀, χ₀ exactly with 0 peaks. Seeds paired across conditions: seed = 20261029 + 1000·(noise index) + replicate.

**Grid.** 200 seeds × 10 cells ((a), (b), (c), (c′), (d) × 3 heights, (e) × 3 heights) × 2 noise levels × 2 modes = 8,000 fits. Success: knee 99.95%, fixed 100%. The 2 knee failures are one spectrum, replicate 163 at nlv 0.02 (seed 20261192), in (a) and (b); with paired seeds those two differ only by a constant offset.

Medians at nlv = 0.02 (nlv = 0.05 in `tables/sim_summary.csv`; same picture, wider spread):

| condition | raw delta | periodic delta, knee | periodic delta, fixed | fixed − knee (paired) |
|---|---|---|---|---|
| (a) baseline | 16.7 | 0.021 | 0.332 | +0.321 |
| (b) offset only | 116.0 | 0.021 | 0.332 | +0.321 |
| (c) exponent only | 14.0 | 0.022 | 0.458 | +0.444 |
| (c′) exponent only, knee freq held | 10.4 | 0.021 | 0.423 | +0.414 |
| (d) peak, true 0.599 | 36.8 | 0.559 | 1.129 | +0.571 |
| (e) all three, true 0.599 | 212.9 | 0.560 | 1.233 | +0.670 |

**What held.**
- Raw delta rises 7× under offset only, as predicted.
- Knee mode: the peak-model periodic delta stays inside its no-peak 95% spread under offset only (94.5%) and exponent only (96%). The quantity does what the paper says it does, when the model matches the generator.
- Knee mode recovers the injected peak with a small downward bias: 0.285 of 0.300, 0.559 of 0.599, 1.114 of 1.198 (−5 to −7%), the same in (d) and (e).

**What did not.**
- *(c) exponent only lowers raw delta (16.7 → 14.0, −16%); the brief predicted a rise.* In b − log10(k + f^χ), f^χ = 1 at 1 Hz for every χ, so with b fixed an exponent increase pivots at 1 Hz and lowers power at every bin above it. "Steeper exponent → more delta" assumes a pivot above the delta band. In the AS-like case the offset dominates (b and c together give +483%). With k fixed, (c) also moves the knee frequency from 3.56 to 2.94 Hz.
- *Fixed mode fails the separation.* Fit to knee-shaped spectra, fixed mode puts ~0.33 of spurious periodic delta into the baseline (the misfit below the knee becomes a delta "peak"; exponent 1.62 vs true 2.21), and an exponent-only change moves it out of its own baseline spread (40% inside at nlv 0.02). It over-recovers the injected peak by 0.2–0.5.

**Primary-estimand preview with known truth.** Fixed minus knee periodic delta is +0.32 with no peak at all, +0.44 when only the exponent changes, +0.57 to +0.67 with the AS-sized peak. Under this generator, the knee-vs-fixed difference is not small, and it is largest exactly where an exponent difference and a delta peak coincide, as in AS.

## 2026-11-03 — (c′) added; reading declared before the rerun

(c′): exponent steepened by Δχ with k = k₀^(χ₁/χ₀), so the knee frequency stays at k₀^(1/χ₀) = 3.56 Hz (in (c) it falls to 2.94 Hz). b stays b₀, so the low-frequency plateau b − log10 k drops by log10(k₁/k₀) = (Δχ/χ₀)·log10 k₀ ≈ 0.21; that is what holding the knee frequency costs in this parameterization.

Reading, written before the run: if the paired fixed − knee periodic delta in (c′) sits near (a)'s +0.32, the knee moving toward the band drives the fixed-mode inflation in (c); if it tracks (c)'s +0.44, the steeper slope does; in between, both contribute.

**Result.** Paired fixed − knee periodic delta, median (IQR):

| | nlv 0.02 | nlv 0.05 |
|---|---|---|
| (a) baseline | +0.321 (0.292–0.348) | +0.308 (0.264–0.353) |
| (c) exponent, knee 3.56 → 2.94 Hz | +0.444 (0.414–0.496) | +0.411 (0.378–0.480) |
| (c′) exponent, knee held at 3.56 Hz | +0.414 (0.381–0.456) | +0.407 (0.360–0.469) |

(c′) tracks (c), not (a): about three quarters of the (c) increase at nlv 0.02 (0.093 of 0.123) and all of it at nlv 0.05. Under the reading declared above, the steeper slope drives the fixed-mode inflation and the knee moving toward the band adds little (~0.03). The IQRs overlap, so "little" is the honest size. PROTOCOL's expectations name both lower knee frequencies and steeper exponents; the simulation puts most of the weight on the exponent.

## 2026-11-03 — Synthetic validation, part 2: parameter recovery, peak detection, noise 0.05

Same run (`tables/sim_summary.csv`, columns `*_err_*`, `frac_any_peak_below_6hz`, `frac_peak_hit_3hz`). Full rerun from scratch with the final code, 12 min 4 s: the 54 rows shared with the earlier table (all of (a)–(e)) are identical in every column, so the earlier hand-filled `true_k` and the summary rebuild changed nothing.

**Knee-mode parameter recovery** (estimated − true, median; nlv 0.02 / 0.05). No-peak cells (a)–(c′): offset −0.014 / −0.035; exponent −0.001 to −0.002 / −0.004 to −0.005; knee frequency −0.02 to −0.03 Hz / −0.05 to −0.06 Hz (IQR about ±0.08 Hz / −0.25 to +0.09 Hz). Peak cells (d)/(e): offset and exponent within ±0.03, knee frequency −0.12 to +0.03 Hz. Knee mode recovers the generator's aperiodic parameters with no bias worth naming. Fixed mode, for reference: exponent −0.59 at baseline (fits 1.62 to a true 2.21), −0.68 in (c′), −0.41 to −0.63 with peaks.

**Peak detection, no-peak cells (a)–(c′): false positives.**

| | any peak < 6 Hz, knee | any peak < 6 Hz, fixed | peak at 3 ± 0.5 Hz, knee | peak at 3 ± 0.5 Hz, fixed |
|---|---|---|---|---|
| nlv 0.02 | 75–77% | 39–81% | 29–30% | 0.5–2.5% |
| nlv 0.05 | 67–68% | 75–86% | 24–25% | 5–12% |

Knee mode fits a small noise peak below 6 Hz in about three quarters of no-peak spectra, and one at 3 ± 0.5 Hz in about a quarter. They are small: the knee-mode periodic delta in those cells has median 0.021 (nlv 0.02) and 0.042–0.050 (nlv 0.05). This is what `peak_threshold` 2.0 with `min_peak_height` 0 and no peak cap does on white noise. The rate should be quoted next to any real-data peak count (secondary (d)).

Fixed mode rarely puts a peak at 3 Hz, yet its baseline periodic delta is 0.33. Spot check, 40 seeds per noise level at baseline (not in the table): the peak contributing most to 2–4 Hz is centered at a median of 5.1–5.2 Hz (5–95%: 3.7–10 Hz), bandwidth ~3.6 Hz, height ~0.37. Fixed mode fills the knee bend with one broad Gaussian near 5 Hz, and its lower tail is the spurious "delta". The 1–6 Hz false-positive rate understates this, because the peak's center is often above 6 Hz.

**Peak detection, (d)/(e): recovery.** A fitted peak at 3 ± 0.5 Hz in 100% of knee fits at both noise levels and all three heights; fixed 96–100% (lowest 96%, A/2 at nlv 0.05).

**Does it hold at nlv 0.05?** Yes, for every claim above.
- Knee mode stays inside its own no-peak 95% spread: (b) 95.0%, (c) 95.0%, (c′) 94.0% (nlv 0.02: 94.5, 96.0, 94.5).
- Knee recovery, raw median vs truth: 0.311 of 0.300, 0.575 of 0.599, 1.125 of 1.198 (+4% to −6%; at 0.02, −5% to −7%). The no-peak baseline is higher at 0.05 (0.050 vs 0.021); recovered minus baseline is −13% to −10% of truth at 0.05 and −12% to −9% at 0.02. Same downward bias at both noise levels.
- Fixed mode, inside its own baseline spread: (c) 72.5%, (c′) 77.5% (nlv 0.02: 40%, 56%).
- Fixed − knee: +0.31 baseline, +0.41 (c), +0.41 (c′), +0.46 (d), +0.58 (e); all lower than at 0.02 by 0–0.11.

One non-monotone result at both noise levels: at 2A the paired gap falls back (d: +0.33 / +0.26; e: +0.48 / +0.38). With a large peak at 3 Hz, the fixed fit uses it in place of the broad 5 Hz shoulder. So the fixed-mode excess is not additive in peak size; it is largest for peaks of the size actually reported for AS (A and A/2).

## 2026-11-04 — Fit pass over the 344 spectra

`python src/fit_spectra.py` (11 s): knee and fixed fits with `fit_one` for every Fig 1a spectrum → `data/features.csv`, `data/peaks.csv`, `data/fits.npz` (all local). Output, in full:

```
converged, knee: TD 185/185, AS 159/159
converged, fixed: TD 185/185, AS 159/159
knee_flag: TD 17/185, AS 8/159
```

`knee_flag` = knee ≤ 0 (no knee frequency) or knee frequency outside 0.5–20 Hz. specparam 2.0.0rc2 fits the knee unbounded (`ap_bounds` ±inf), so "knee at bound" (gate item 4) cannot occur; the flag covers the implausible-knee half. Canary: TD knee medians b₀ = 2.359707, k₀ = 16.518708, χ₀ = 2.208844, and 9 knee ≤ 0 + 8 out of range, identical to the synthetic baseline entry, so `fit_one` is unchanged by the new `return_components` option.

## 2026-11-04 — QC overlay selection and age bins, declared before running `src/qc.py`

Ten overlays per group, chosen by rule, not by browsing:

1. **Highest raw delta**: the recording with the largest `raw_delta` (our band rule on Fig 1a, identical to the reported delta to 4.7e-7).
2. **Knee-flagged**: among the group's `knee_flag` recordings, the one with the lowest knee-mode R². If the group had none, the recording whose knee frequency is farthest from the group median on a log scale.
3. **Named**: AS154 (sharp 2.5 Hz peak) for AS; TD133 (footnoted offset) for TD.
4. **Collisions**: if a rule picks a recording already chosen, it takes the next one in the same ordering. Order of filling: named, then highest raw delta, then knee-flagged.
5. **Seven random**: `np.random.RandomState(20261029)`, one generator per group (TD first, then AS, each freshly seeded), `choice(remaining, 7, replace=False)` over the group's remaining labels in Fig 1 order.

Age bins for the convergence and knee-flag counts, left-closed: [0, 3), [3, 6), [6, 9), [9, 12), [12, 16) years (n TD/AS: 48/37, 54/48, 48/39, 23/23, 12/12).

QC does not compute or plot any summary of fixed-mode periodic delta or the fixed − knee difference across recordings; the overlays show both fits per recording, as gate item 5 requires.

## 2026-11-04 — QC outputs (gate items 2–5)

`python src/qc.py` → `tables/qc_fit_quality.csv`, `tables/qc_counts.csv`; local: `data/qc/knee_flagged.csv`, `data/qc/hist_r2_error.png`, `data/qc/overlay_selection.csv`, `data/qc/overlay_<label>.png` (20), `data/qc/overlays_TD.png`, `data/qc/overlays_AS.png`.

**Convergence** 344/344 in both modes, every age bin.

**Knee flags** (all are knee ≤ 0 or knee frequency < 0.5 Hz; none above 20 Hz):

| age bin | TD flagged / n | AS flagged / n |
|---|---|---|
| [0, 3) | 13 / 48 | 3 / 37 |
| [3, 6) | 3 / 54 | 4 / 48 |
| [6, 9) | 0 / 48 | 1 / 39 |
| [9, 12) | 1 / 23 | 0 / 23 |
| [12, 16) | 0 / 12 | 0 / 12 |
| all | 17 / 185 (9 knee ≤ 0) | 8 / 159 (6 knee ≤ 0) |

The TD flags are concentrated under age 3 (13 of 17; 11 of them under 1.1 years). These are spectra with no bend for the knee to fit, so the knee fit collapses toward a straight line.

**R² and error** (median [5th percentile, min]):
- Knee: TD 0.9970 [0.990, 0.951], AS 0.9970 [0.983, 0.921].
- Fixed: TD 0.9928 [0.983, 0.953], AS 0.9942 [0.981, 0.920].
- Error medians: knee TD 0.029, AS 0.040; fixed TD 0.041, AS 0.046.
- The two lowest R² in each mode are TD150 and AS151, both knee-flagged.

**Overlays, one sentence each** (selection per the rule declared above; `data/qc/overlay_selection.csv`).

TD:
- TD133 (named): knee ≤ 0, so both aperiodic fits are the same straight line; both peak models follow a large 2–8 Hz bump and a 20–50 Hz residual. Sane given a spectrum with no bend.
- TD23 (highest raw delta): the knee fit follows the low-frequency plateau and is sane. The fixed line overshoots at 1 Hz and undershoots 2–8 Hz, so its peak model adds a broad 2–8 Hz shoulder under the alpha peak.
- TD150 (knee-flagged, lowest knee R² 0.951): not sane in the delta band. Both fits draw a straight line through a 1.5–6 Hz bump (raw residual up to +0.28), and neither peak model fits it; their only peak is at 30 Hz. See the edge-rule note below.
- TD93: knee frequency 0.85 Hz, so the two fits are nearly identical; peak models agree. Sane.
- TD36: knee sane (plateau below 4 Hz). Fixed overshoots at 1 Hz and its peak model carries a broad 3–8 Hz shoulder and inflated alpha and beta.
- TD167: knee sane; fixed similar, with a shoulder from 4 Hz.
- TD182: knee sane; fixed adds a broad 3–7 Hz shoulder and a taller alpha peak.
- TD112: both reasonable; fixed has a slightly larger 2–6 Hz bump; the knee peak model fits small ripples at 30–50 Hz.
- TD156: both sane and close.
- TD59: knee sane; fixed adds a 2–6 Hz shoulder.

AS:
- AS154 (named): the knee peak model captures the 2.5 Hz peak but under-fits its height (about 1.25 vs residual about 1.55), as noted on 2026-11-01. The fixed peak model is taller and adds a 4–20 Hz plateau of peaks. Knee sane apart from the known under-fit.
- AS105 (highest raw delta): sane in shape, but the knee aperiodic plateau (knee frequency 3.59 Hz) runs through the 2–4 Hz excess, so part of the delta bump is assigned to the aperiodic component. Fixed assigns all of it to peaks.
- AS151 (knee-flagged, lowest knee R² 0.921): not sane in the delta band. Both fits are straight lines, and neither peak model fits any peak, despite a 1.5–4 Hz bump of +0.85 in the residual. Periodic delta is 0 in both modes.
- AS58: both sane; fixed larger at 2–5 Hz.
- AS116: both sane; fixed larger at 2–5 Hz.
- AS146: knee sane; fixed about twice as high in the delta band.
- AS20: knee sane with a small delta peak model; fixed puts a broad 3–8 Hz peak there.
- AS45: knee sane, with its plateau ending near 2 Hz; fixed about twice as high from 2 to 8 Hz.
- AS41: both sane.
- AS22: knee sane; fixed adds a broad shoulder from 2 Hz.

**Edge rule (diagnostic, not a change).** specparam drops a candidate peak whose center lies within 1 × its bandwidth of the edge of the fit range (`_bw_std_edge = 1.0`, the default the paper also ran). Refitting TD150 and AS151 in knee mode with `_bw_std_edge = 0` keeps a peak at 4.16 Hz (TD150) and 2.86 Hz (AS151), and the peak-model delta goes from 0.000 to 1.02 and 1.76. So a broad low-frequency bump next to the 1 Hz edge can be reported as zero periodic delta under the paper's settings. The analysis keeps the default. How often this happens across the 344 is a candidate for secondary (d); not counted here.

**Gate items, from these outputs** (verdict in the gate entry below):
- Item 2 as written ("both fits succeed for ≥ 95%, no systematic failure by age or group"): convergence passes (100%). Under PROTOCOL's Inclusion definition of a failed fit (no convergence or knee frequency outside 0.5–20 Hz), success is 319/344 = 92.7%: TD 168/185 = 90.8%, AS 151/159 = 95.0%. TD failures are concentrated under age 3.
- Item 3: reported above; no threshold applied.
- Item 4: counted (table above) and listed (`data/qc/knee_flagged.csv`).
- Item 5: 20 overlays saved and inspected.

## 2026-11-04 — Validation table (secondary a): inclusion declared before running `src/validation.py`

Knee refit vs reported (Fig 1b,d,f,g), per group: exponent, offset, periodic delta (refit `periodic_delta_peak_knee` vs reported periodic delta). Difference = refit − reported. All 344 recordings included (all converged); knee-flagged recordings stay in. A one-line sensitivity excludes the 25 knee-flagged. Measures: n, Pearson r, mean bias, SD of the difference, 95% limits of agreement (bias ± 1.96 SD), and the median difference.

Prediction on record (from the Oct 25 check): TD exponent within ~0.05–0.10, and AS about 0.15 off, "high". The one recorded data point is TD1: refit 1.935 vs reported 1.981, so refit − reported = −0.046. "AS high" is read here as the AS gap being about 0.15 in magnitude. The sign will be reported as found.

## 2026-11-04 — Validation table (secondary a): result

`python src/validation.py` → `tables/validation.csv`, `tables/deviations.csv`. Knee refit − reported, analysis set (all 344; sensitivity without the 25 knee-flagged in brackets):

| | n | r | mean bias | 95% LoA | median \|d\| | \|d\| > 0.15 |
|---|---|---|---|---|---|---|
| TD exponent | 185 | 0.986 | −0.012 [−0.004] | −0.22 to +0.19 | 0.063 | 13% |
| AS exponent | 159 | 0.929 | −0.002 [+0.006] | −0.35 to +0.35 | 0.079 | 20% |
| TD offset | 185 | 0.985 | +0.040 [+0.048] | −0.24 to +0.32 | 0.076 | 28% |
| AS offset | 159 | 0.932 | +0.055 [+0.062] | −0.46 to +0.57 | 0.115 | 38% |
| TD periodic delta | 185 | 0.819 | −0.017 [−0.009] | −0.25 to +0.22 | 0.051 | 9% |
| AS periodic delta | 159 | 0.914 | −0.125 [−0.116] | −0.51 to +0.26 | 0.114 | 36% |

**Against the prediction.**
- TD exponent "within ~0.05–0.10": held as a typical size (median |d| 0.063; mean bias −0.012). The per-recording spread is wider than the prediction implies (LoA about ±0.2), and that spread is the finding. TD1 (−0.046) sits inside the TD interquartile range (−0.073 to +0.056).
- AS exponent "about +0.15 high": not held as a direction. Mean −0.002, median +0.015, refit above reported in 55% and below in 45%. Read as a size, partly: median |d| 0.079, and 20% of AS recordings differ by more than 0.15.
- So the exponent pattern that would have become the email's second question did not appear across all 159.

**Not predicted: AS periodic delta.** The knee refit's periodic delta is lower than reported in AS: mean −0.125, median −0.104, 36% of recordings off by more than 0.15, and the same without the knee-flagged (−0.116). In TD it is essentially zero (−0.017). This was not predicted. It is the stronger candidate for the email's second question, in place of the exponent gap. (TD r = 0.82 reflects TD's small range of periodic delta, mean 0.13; compare the bias and LoA, not r.)

**Not investigated here.** Same specparam version and settings, and the paper's own per-recording aperiodic fit (Fig 1e) and peak model (Fig 1c) are in the supplement. Comparing our knee aperiodic fit with Fig 1e recording by recording would show whether the fitting or the inputs differ. That is outside this phase; it is the open question the email would ask.

**Offsets.** Refit about +0.04 (TD) and +0.055 (AS) above reported, same sign in both groups. Recorded in `tables/deviations.csv`.

**What has now been seen.** The PROTOCOL disclosure (2026-11-03) stays as written; since then, knee-mode aggregates per group for exponent, offset and periodic delta have been computed under secondary (a), and the 20 QC overlays show both modes per recording. Still unseen: any aggregate of fixed-mode periodic delta, and any fixed − knee difference across recordings or groups.

**2026-11-05 — Prediction resolved.** The Oct 25 prediction (“AS about +0.15 high”) rested on three recordings and is retracted as a direction; as a size it partly holds (median |d| 0.079, 20% beyond 0.15). The finding that replaces it: no mean bias in either group but per-recording limits of agreement of about ±0.2 (TD) and ±0.35 (AS), i.e. knee-mode parameters are not identified to better than that on a single spectrum. The AS periodic-delta gap (−0.125) is the open question and is the subject of the Fig 1e diagnostic, with readings declared before the run:

*Fig 1e diagnostic (PROTOCOL secondary (a), continued), readings declared before the run.* This replaces a draft written earlier the same day (cuts 1e-5 / 1e-4), which was never run.

Fig 1e is the paper's aperiodic component per recording, in linear units. If it came from one specparam knee fit on the region-averaged spectrum, log10(Fig 1e) is exactly b − log10(k + f^χ), apart from spreadsheet rounding (7 significant figures, about 2e-7 in log10). Per recording: `scipy.optimize.least_squares` fit of that curve to log10(Fig 1e) over the 98 bins, starting from the reported offset and exponent with k = 3^χ (a 3 Hz knee); then the largest absolute residual, and fitted (b, χ) vs the reported offset (Fig 1g) and exponent (Fig 1f).

| outcome | max residual | fitted (b, χ) vs reported | meaning |
|---|---|---|---|
| (A) | < 1e-4 | equal within 1e-3 | Fig 1e is a single fit, and the reported parameters are that fit. The refit's ±0.2 scatter then comes from the inputs or the optimizer, not from aggregation. Next: check whether the refit's aperiodic curve differs from Fig 1e mainly in 2–4 Hz, and compare the refit's peak model with Fig 1c across the full spectrum. |
| (B) | > 1e-3, structured | n/a | Fig 1e is not one knee curve: an average of regional fits (or smoothed). An average of three knee curves is not a knee curve. This explains scatter without bias in exponent and the −0.125 AS periodic-delta gap (averaging spectra blurs peaks; a fit on the average finds less peak power than the average of per-region peak fits). The email's second question becomes a request to confirm. |
| (C) | < 1e-4 | differ | Fig 1e and Fig 1f/1g come from different computations; ask which the reported features are. |

Operational details fixed now:
- "Structured" means the residual is smooth across frequency, not rounding noise: lag-1 autocorrelation of the residual over the 98 bins ≥ 0.5.
- A recording with max residual between 1e-4 and 1e-3, or above 1e-3 but not structured, is "unclassified" and is listed, not assigned.
- Counts per group (TD, AS) in each outcome are reported. A mixture is reported as a mixture, not forced into one outcome.

**Two side checks, declared with the readings.**
1. *Edge-rule cases.* Do TD150 and AS151 have nonzero reported periodic delta? If yes, the paper's fit saw the low-frequency bump that the edge rule drops in our refit. If they are 0, the paper's fit dropped it too.
2. *Does the AS gap scale with reported periodic delta?* Blurring predicts it does: refit − reported periodic delta becomes more negative as reported periodic delta grows. Measures, per group: OLS slope of (refit − reported) on reported, with its 95% CI, and Spearman ρ. Because the difference contains −reported, regressing it on reported is biased toward a negative slope by measurement noise alone (regression to the mean). So the same slope is also reported against the pair mean, (refit + reported)/2, as in a Bland–Altman plot. Blurring is supported only if the slope is negative on both axes in AS. A negative slope on the reported axis alone does not count. TD is reported alongside for reference.

## 2026-11-04 — QC gate verdict, item by item

`requirements.txt` (exact pins from the project venv's freeze, minus fooof, which only the smoke tests use) and `scripts/run_all.sh` (supplement sha256 check → `validate.py` → `fit_spectra.py` → `qc.py` → `validation.py` → pytest; `--with-sim` adds the 12-minute synthetic grid).

1. **Loader assertions pass: PASS.** pytest 19/19 green (13 loader, 2 config, 4 simulate). `load_all`'s assertions also run inside every fit pass.
2. **Both fits succeed for ≥ 95%, no systematic failure by age or group: FAILS as written.**
   - Convergence is 344/344 in both modes, in every age bin.
   - PROTOCOL's Inclusion section defines a failed fit as no convergence *or* knee frequency outside 0.5–20 Hz. By that definition, knee mode succeeds for 319/344 = 92.7% (TD 168/185 = 90.8%, AS 151/159 = 95.0%).
   - The failures are systematic by age in TD: 13 of the 17 TD flags are under age 3 (13/48 = 27% of that bin, vs 3/37 in AS).
   - Consequence if the criterion stands: 25 recordings leave the paired comparison, 13 of them TD under 3.
   - Not resolved here;
3. **R² distributions reported, no borrowed threshold: PASS.** `tables/qc_fit_quality.csv`; histograms in `data/qc/hist_r2_error.png`.
4. **Knee-at-bound or implausible knee frequency counted and listed: PASS.** Knee-at-bound cannot occur (knee unbounded in specparam 2.0.0rc2); implausible knee: 25, counted by group and age bin in `tables/qc_counts.csv`, listed in `data/qc/knee_flagged.csv`.
5. **Ten overlays per group inspected and saved: PASS**, with two findings.
   - Selection was by the rule declared before running. The 20 PNGs are in `data/qc/`, with one sentence each in the QC entry.
   - 18/20 fits look sane.
   - TD150 and AS151 (both knee-flagged) are not sane in the delta band: a low-frequency bump is dropped by specparam's edge rule and reported as zero periodic delta.
6. **Whole thing reruns from a clean environment: PASS.**
   - Clean copy of `src/`, `tests/`, `scripts/`, `refs/`, `config.yaml`, `pyproject.toml`, `requirements.txt` (no `data/`, `tables/` or venv). Fresh `venv` from Homebrew Python 3.13.7, `pip install -r requirements.txt` (freeze matched the file exactly), then `bash scripts/run_all.sh`: exit 0.
   - Every output is byte-identical to the project copy: `features.csv`, `peaks.csv`, `validate_delta.csv`, all QC CSVs and the 23 QC PNGs, and the five tables written there (`delta_rules`, `qc_counts`, `qc_fit_quality`, `validation`, `deviations`). All 9 arrays in `fits.npz` are equal.
   - Scope: same machine and base interpreter; the synthetic grid (`--with-sim`) was not part of this rerun. Its own from-scratch reproduction is in the 2026-11-03 entry.

**Verdict:** 5 of 6 pass; item 2 fails as written. The primary analysis is not run until item 2 is decided.

**2026-11-05 — Item 2 resolved.** Decision: item 2 means convergence, as v2.2’s Inclusion section already states; the item’s wording is corrected in v2.3. Convergence is 344/344 in both modes and every age bin: PASS. The 25 knee-flagged recordings stay in the primary analysis and are dropped only in the predeclared sensitivity, because excluding them conditions on the knee-mode outcome, one arm of the comparison, and would push the paired median away from zero on exactly the recordings where the two models coincide. Gate: 6 of 6. The primary analysis may run after the Fig 1e diagnostic above.

## 2026-11-05 — Fig 1e diagnostic and side checks: results

`python src/fig1e_diagnostic.py` → `tables/fig1e_diagnostic.csv`, `tables/periodic_gap_scaling.csv`; per recording in `data/fig1e_fits.csv` (local). Added to `scripts/run_all.sh` after `validation.py`.

**Fitter check first** (`tests/test_fig1e.py`). The test rounds to 7 significant figures; the supplement is rounded to 6 decimal places, coarser at high frequencies (4 significant figures at 49.5 Hz). The worst rounding bound in log10 is 3e-5, below every max residual observed, so the test’s gentler rounding does not affect the outcome.

**Outcomes, by the readings declared above.** All 344 fits converged.

| | (A) | (B) | (C) | unclassified |
|---|---|---|---|---|
| TD (185) | 0 | 120 | 0 | 65 |
| AS (159) | 0 | 91 | 0 | 68 |

- **Precision correction.** The readings above say Fig 1e is stored to 7 significant figures. That holds at the low end (17.52503 at 1 Hz); at the high end it is 6 decimal places (0.007249, the smallest value, at 49.5 Hz). So rounding error in log10 runs from about 2e-7 at 1 Hz to about 3e-5 at the smallest values. Both are below the 1e-4 cut, so the declared cuts are unaffected.
- **No recording is a single knee curve.** The smallest max residual of the 344 is 2.0e-4 (AS109), at 1.0 Hz, where rounding is about 2e-7, so about 1,000 times rounding. The largest residual is at 1.0 or 1.5 Hz in 341 of 344 recordings; the median is 7.4e-3 (5th–95th percentile 1.0e-3 to 3.3e-2). (A) and (C) are rejected for every recording, and that holds whatever "structured" means.
- **211 are (B)** by the declared criterion (max residual > 1e-3 and lag-1 autocorrelation ≥ 0.5).
- **133 are unclassified.** 17 have a max residual between 1e-4 and 1e-3. 116 exceed 1e-3 but have autocorrelation < 0.5.
- **Where the misfit sits (observation, not a reclassification).** The largest residual is at 1.0 Hz in 310 of 344 recordings, and at 1.5 Hz in 31 more. In the unclassified ones it is concentrated in the lowest three or four bins (1–2.5 Hz), often with a sign change between 1.0 and 1.5 Hz (AS33: −0.020, +0.022, +0.014 at 1.0, 1.5, 2.0 Hz), and near zero above. A difference concentrated at the low edge does not register in a lag-1 autocorrelation over 98 evenly spaced bins, so the declared "structured" test undercounts it. The declared counts stand as reported. Whether to reclassify with a different structure measure is a separate, labeled-as-post-hoc decision.
- **Fitted vs reported parameters.** Median |b − offset| 0.0035, median |χ − exponent| 0.0023, largest 0.23 and 0.15 (AS95). The best single curve through Fig 1e is close to the reported parameters but not equal to them.

**Reading.** Fig 1e is not the curve of one knee fit, for any recording. That rules out (A) and (C) and is consistent with (B): an average of fitted curves, or a smoothed curve. The data do not say which kind of averaging. The misfit is largest at 1–2.5 Hz, where knee curves with different knees differ most.

**Side check 1: edge-rule cases.** Reported periodic delta is TD150 0.795 and AS151 1.448; our knee refit gives 0 for both. The paper's fit saw the low-frequency bumps that the edge rule drops in a single fit of Fig 1a. That is consistent with fits made on spectra other than Fig 1a (for example per region, where the bump need not sit at the 1 Hz edge), but does not show it.

**Side check 2: does the gap scale with reported periodic delta?** Refit − reported, OLS slope [95% CI], Spearman ρ:

| | vs reported | vs pair mean |
|---|---|---|
| AS | −0.075 [−0.139, −0.010], ρ −0.31 | +0.013 [−0.054, +0.080], ρ −0.17 |
| TD | −0.254 [−0.330, −0.178], ρ −0.58 | −0.102 [−0.194, −0.010], ρ −0.47 |

- **AS: blurring not supported by the declared rule.** The slope is negative against reported only, which regression to the mean produces by itself; against the pair mean it is flat. The AS gap of −0.125 is an offset that does not grow with reported periodic delta.
- **TD:** the slope is negative on both axes. This was not predicted, and the mean TD gap is about 0 (−0.017).

**Where this leaves the open question.** Fig 1e is not a single knee fit, so the reported features come from a computation other than one knee fit of the Fig 1a spectrum. That is the email's second question, now a request to confirm: were the aperiodic component and the reported features computed per region (or per channel) and then averaged? The AS periodic-delta gap is not explained by the blurring mechanism as tested; it stays open.

**Gate item 6, extended.** The clean-environment rerun (same fresh venv, clean copy updated with `src/fig1e_diagnostic.py`, `tests/test_fig1e.py` and the new `run_all.sh`) exits 0 with 24 tests passing. `fig1e_diagnostic.csv`, `periodic_gap_scaling.csv` and `data/fig1e_fits.csv` are byte-identical to the project copy, and so are all earlier outputs.

Per the gate resolution, the primary analysis may now run.

## 2026-11-14 — Gap subsets and Fig 1e cross-tabs (post hoc, requested after the Nov 5 results)

`python src/gap_subsets.py` → `tables/periodic_gap_subsets.csv`, `tables/fig1e_crosstabs.csv`. Neither test was declared before the Nov 5 results; both are labeled post hoc.

**(i) Gap (knee refit − reported periodic delta) by refit peak-model delta.**

| | n | mean gap | median gap [IQR] | mean reported | share of the group's summed gap |
|---|---|---|---|---|---|
| AS, refit = 0 | 4 | −0.936 | −0.812 [−1.046, −0.702] | 0.936 | 19% |
| AS, refit > 0 | 155 | −0.104 | −0.098 [−0.171, −0.035] | 0.701 | 81% |
| TD, refit = 0 | 0 | — | — | — | — |
| TD, refit > 0 | 185 | −0.017 | −0.006 [−0.055, +0.045] | 0.133 | 100% |

- The exact-zero split tests nothing in TD: no TD refit is exactly 0. A fitted peak anywhere leaves a small Gaussian tail in 2–4 Hz, so "> 0" includes near-zero values (89 TD and 23 AS recordings are between 0 and 0.01).
- **Added alongside (my split, not requested): < 0.01 vs ≥ 0.01.** AS < 0.01: n 27, mean −0.287 (39% of the AS gap); AS ≥ 0.01: n 132, mean −0.091 (61%). TD < 0.01: n 89, mean −0.028, median +0.009; TD ≥ 0.01: n 96, mean −0.007.
- **Reading.** The four AS recordings with refit exactly 0 (AS68, AS119, AS150, AS151) carry 19% of the AS gap. All four have a knee-mode raw-residual periodic delta above 0.3 (0.61–1.11), so all four count as edge-rule drops under PROTOCOL secondary (d). Knee-mode count of that measure: AS 4, TD 0; the fixed-mode count waits for the primary analysis. Excluding them leaves an AS gap of −0.104, and among recordings with a refit delta ≥ 0.01 it is −0.091, against −0.007 in TD. Edge-rule drops add to the AS gap but do not account for most of it. Most of it sits in recordings where the refit does fit delta peaks and finds less than reported.

**(ii) Fig 1e status vs refit knee flag and vs age < 3, 2 × 2 per group** (Fisher exact tests, descriptive):

| group | status | vs | status & yes | status & no | not status & yes | not status & no | odds ratio | p |
|---|---|---|---|---|---|---|---|---|
| TD | unclassified | knee flag | 15 | 50 | 2 | 118 | 17.7 | 3e-6 |
| TD | unclassified | age < 3 | 41 | 24 | 7 | 113 | 27.6 | 4e-17 |
| TD | flip | knee flag | 16 | 45 | 1 | 123 | 43.7 | 6e-8 |
| TD | flip | age < 3 | 44 | 17 | 4 | 120 | 77.6 | 8e-24 |
| AS | unclassified | knee flag | 7 | 61 | 1 | 90 | 10.3 | 0.02 |
| AS | unclassified | age < 3 | 25 | 43 | 12 | 79 | 3.8 | 6e-4 |
| AS | flip | knee flag | 8 | 58 | 0 | 93 | ∞ | 7e-4 |
| AS | flip | age < 3 | 27 | 39 | 10 | 83 | 5.7 | 2e-5 |

"Unclassified" is the declared Fig 1e outcome against (B) (there are no (A) or (C)). "Flip" means the single-curve residual changes sign between 1.0 and 1.5 Hz.

- **Reading.** The edge-localized misfit belongs to the youngest children and to recordings whose refit knee is flagged. In TD, 44 of the 48 recordings under 3 flip, against 17 of 137 older; 16 of the 17 knee-flagged TD recordings flip. AS shows the same direction, more weakly.
- These are the spectra with little or no bend above 1 Hz. If Fig 1e averages several knee curves (outcome B), curves whose knees sit at or below 1 Hz would differ from a single knee curve mainly in the lowest bins. That is the pattern seen, and it is the pattern the lag-1 structure test misses. This is consistent with (B); it does not establish it.
- Both status variables and the knee flag are partly functions of the same spectral shape, so the associations are not independent evidence.

## 2026-11-14 — Gate item 6 with the synthetic grid: 6/6, synthetic reproduced, 2026-11-14, commit 7a2bf66

- **What ran.** `git archive 7a2bf66` (clean tree; HEAD at the time) exported to an empty directory; a fresh `venv` from Homebrew Python 3.13.7 at `./.venv` inside it (the path `run_all.sh` now uses); `python -m pip install -r requirements.txt` (freeze identical to the file); `bash scripts/run_all.sh --with-sim`.
- **Result.** Exit 0 in 12 min 46 s. Supplement sha256 OK; 344/344 converged in both modes; knee flags TD 17, AS 8; Fig 1e outcomes as on Nov 5; 8,000 synthetic fits; 24/24 tests.
- **Committed outputs.** All 8 tables in `tables/` and `figures/fig_synthetic_validation.png` regenerated byte-identical to their committed versions at 7a2bf66, including `tables/sim_summary.csv` and the figure.
- **Local outputs.** Every file the run wrote to `data/` is byte-identical to the project's local copy: `features.csv`, `peaks.csv`, `fig1e_fits.csv`, `sim_spectra.csv`, `validate_delta.csv`, the QC CSVs and all 23 QC PNGs. All 9 arrays in `fits.npz` are equal.
- **Gate: 6 of 6**, with item 6 now covering the synthetic grid as well as the analysis pipeline.
- **Not covered by the hash.** `src/gap_subsets.py` and its two tables (2026-11-14, uncommitted when this ran) are not in 7a2bf66. Run in the same fresh environment afterwards, they reproduce the project's tables byte-identically, but that check is not pinned to a commit.

## 2026-11-14 — Primary analysis and secondaries (a)–(e)

`python src/analysis.py` (6 s; in `run_all.sh` after `qc.py`) → `tables/primary.csv`, `tables/secondary_a.csv` … `secondary_e.csv`, `tables/group_contrast_by_mode.csv`; `figures/fig_primary_paired.png`, `figures/fig_knee_flags_age.png`, `figures/fig_group_contrast_by_mode.png`. Local only: `data/fig_primary_paired_points.png`. The label of every row (predeclared / exploratory / post hoc) was fixed in the script's docstring before the first run.

**Spec vs PROTOCOL, and other choices.**
- *Inputs.* `features.csv` is wide (one row per recording, `_knee` / `_fixed` columns), not one row per recording × mode; used as is. (a) and the gap decomposition in (d) compare with the paper's reported values, which are not in `features.csv`; those two read them through `load_supplement.load_reported`. Nothing else is read, nothing is refit.
- *Cluster bootstrap.* Participant IDs are in `features.csv` (Fig 2 ID: 185 TD participants, one recording each; 95 AS participants), so the predeclared sensitivity runs. Cluster key = group + ID. TD's cluster CI equals its recording CI by construction (185 singletons); it differs only for AS and overall.
- *(b).* The spec's Δexponent and Δoffset (fixed − knee) by group are in `secondary_b.csv`, labelled exploratory. PROTOCOL (b) is the AS-vs-TD Welch Cohen's d per mode for periodic delta and exponent; it is in `group_contrast_by_mode.csv`, labelled predeclared. `secondary_b.csv` also has knee-flagged and unflagged subsets, added to test the 2026-11-05 expectation.
- *(c).* Slopes on log2(age), i.e. per doubling of age (PROTOCOL says log age; the spec says age). The slope of the paired difference equals the fixed − knee slope difference, so the Δexponent slope is PROTOCOL's "slope difference by paired bootstrap". PROTOCOL (c)'s per-mode TD slopes for exponent, offset and periodic delta are included. AS rows are exploratory. No GAM.
- *(d).* PROTOCOL's convergence, fit-error and low-frequency peak counts are included alongside the spec's three blocks. The gap decomposition is post hoc (2026-11-14) and labelled so.
- *Group contrast.* The AS − TD median difference per mode, and its change from knee to fixed, are exploratory: PROTOCOL names the contrast only as Cohen's d.
- *Figure.* "d against the pair mean" as a scatter would put per-recording values in a committed figure. The committed version bins the pair mean into quintiles per group (37 TD / 31–32 AS per bin), showing the bin median and IQR of d plus the group median line and CI band. The point scatter is local.
- *Bootstrap.* 10,000 resamples, seed 20261029, percentile 95% CI; each interval uses its own freshly seeded generator, so each CI is reproducible on its own.
- *Cross-checks.* `secondary_a.csv` equals `validation.csv` on all 12 shared rows; the knee-flag counts in `secondary_d.csv` equal `qc_counts.csv`.

**Primary: d = fixed − knee periodic delta (peak model).**

| | n | median d [95% CI] | participant-cluster CI | d > 0 | mean d [CI], descriptive |
|---|---|---|---|---|---|
| all | 344 | +0.197 [0.161, 0.229] | [0.158, 0.233] | 92% | +0.268 [0.241, 0.297] |
| TD | 185 | +0.110 [0.094, 0.129] | same (singletons) | 89% | +0.131 [0.113, 0.150] |
| AS | 159 | +0.426 [0.349, 0.481] | [0.331, 0.482] | 96% | +0.428 [0.384, 0.474] |

Without the 25 knee-flagged recordings (predeclared sensitivity): all +0.217 [0.183, 0.241], TD +0.123 [0.102, 0.144], AS +0.446 [0.394, 0.502]. Excluding them moves the medians up slightly, as expected for recordings where the two fits coincide; the conclusions do not change.

**(e) Same comparison on the raw-residual periodic delta:** all +0.282 [0.254, 0.309], TD +0.205 [0.182, 0.231], AS +0.490 [0.431, 0.553]. Same direction and ordering, larger in TD.

**What the study set out to say (PROTOCOL "can say").**
- *Is the AS–TD contrast stable under the change of model?* No: it grows. Periodic delta Cohen's d (PROTOCOL (b)) goes from 1.27 [1.08, 1.49] under knee to 1.72 [1.50, 1.98] under fixed; exponent from 0.87 [0.66, 1.11] to 1.53 [1.32, 1.76]. AS − TD median periodic delta goes from +0.472 [0.376, 0.649] to +0.715 [0.619, 0.893], a change of +0.243 [0.135, 0.339] (exploratory).
- *Is the TD age association stable?* No, not for exponent. The TD exponent slope (PROTOCOL (c)) is +0.320 per doubling of age [0.270, 0.370] under knee and +0.078 [0.050, 0.107] under fixed, about three quarters smaller. Offset: +0.447 vs +0.061. Periodic delta is closer: −0.111 vs −0.088, with a slope difference of +0.024 [0.008, 0.039].
- *Exploratory (AS):* the AS exponent age slope changes sign between modes: knee +0.070 [−0.004, +0.149], fixed −0.109 [−0.146, −0.070].

**(b) Fixed − knee aperiodic parameters (exploratory).** Δexponent TD −0.504 [−0.631, −0.439], AS −0.603 [−0.696, −0.531]; Δoffset TD −0.788, AS −0.915. Positive in only 15 recordings, the same 15 for both parameters, 13 of them with knee ≤ 0.

**(d)** Edge-rule count (peak-model delta 0, raw-residual delta > 0.3): knee TD 0, AS 4; fixed TD 0, AS 3. Knee-flag rates by age bin as in the QC entry (`fig_knee_flags_age.png`). Convergence 344/344 in both modes. Mean number of peaks below 6 Hz: knee TD 0.68, AS 1.40; fixed TD 0.78, AS 1.28.

**Expectations, clause by clause.**

PROTOCOL, 2026-11-03 (from synthetic validation):

| clause | result | verdict |
|---|---|---|
| fixed exceeds knee by roughly 0.3 in TD | +0.110 [0.094, 0.129]; 89% of TD recordings d > 0 | direction held; **size not held** (the CI excludes 0.3; about a third of it) |
| … and by more in AS | +0.426 [0.349, 0.481] | held |
| the AS–TD contrast under fixed mode does not vanish and may increase | Cohen's d 1.27 → 1.72; median difference +0.47 → +0.72, change +0.24 [0.14, 0.34] | held; it increased |
| because AS knee frequencies are lower and exponents steeper, for a reason unrelated to oscillations | — | **not tested** by these numbers: they measure the change in the contrast, not its cause |

PROTOCOL, 2026-11-05 (age pattern of knee flags):

| clause | result | verdict |
|---|---|---|
| Δexponent ≈ 0 for knee-flagged recordings | TD +0.028 [−0.008, +0.067] (n 17); AS −0.000 [−0.039, +0.050] (n 8) | held |
| an age trend in Δexponent | TD −0.241 per doubling of age [−0.274, −0.209] | present |
| the flagged recordings may themselves produce that trend | without them: −0.227 [−0.266, −0.191], 94% of the full-sample slope | **not supported**: the trend is in the unflagged recordings |

So two of the four 2026-11-03 clauses held, one failed on size, and one is untested. Of the 2026-11-05 clauses, the first held, the trend is present, and the proposed cause of the trend is not supported.

**Open: the synthetic generator overpredicted the TD gap about 3×.** Its TD-like no-peak baseline gave fixed − knee +0.32; the real TD median is +0.11. The generator has one knee curve, white log10 noise and no real peaks. Not explained here.

**What has now been seen.** As of this entry, fixed-mode aggregates and the fixed − knee contrast have been computed. Nothing predeclared remains unseen.

**Reproduction.** A fresh copy of the tracked and new files (no `data/`) with a new venv from `requirements.txt` and `bash scripts/run_all.sh`: exit 0, 27/27 tests. All 51 outputs are byte-identical to the project's (17 tables, 4 figures, 30 local `data/` files). `sim_summary.csv` and the synthetic figure were copies, not regenerated (no `--with-sim`); their reproduction is the 7a2bf66 entry. This rerun is of uncommitted files, so it is not pinned to a commit.

## 2026-11-16 — Why the synthetic baseline overpredicted the TD gap: design and readings, declared before the run

Post hoc and exploratory. Written after the primary results (2026-11-14) and after printing the TD knee-fit quantiles and bin counts below, before computing any fixed − knee value by bin or fitting any twin.

**Gap to explain.** Synthetic (a), nlv 0.02: paired fixed − knee periodic delta, median +0.321 (200 replicates). Real TD: +0.110. Difference 0.21.

**What differs.** (a) is one knee curve at the TD medians (b₀ 2.36, k₀ 16.5, χ₀ 2.21; knee 3.56 Hz) with no peaks and white log10 noise (SD 0.02). Real TD knee fits spread widely: knee frequency 1.0–6.4 Hz (10th–90th percentile; 9 more have knee ≤ 0 and no knee frequency), χ 1.5–3.0. Their knee-model residual error has median 0.029, and they have 4 fitted peaks (median). Knee fits of the generator's spectra find 11 peaks (median), all noise.

**Step 1, no new fits (`features.csv`).** TD median d by knee-frequency bin, with n per bin:

| bin | < 2 Hz | 2–3 Hz | 3–4.2 Hz | 4.2–6 Hz | ≥ 6 Hz | no knee frequency |
|---|---|---|---|---|---|---|
| n | 43 | 25 | 33 | 47 | 28 | 9 |

Also computed in recordings near the baseline (knee 3.0–4.2 Hz and χ 2.0–2.4; n 13). Each row also gets the median fixed and median knee periodic delta, because in TD the difference of medians (0.228 − 0.015) is twice the median difference.

Reading:
- **Near-baseline median d above 0.215** (closer to the generator's 0.32 than to 0.11): recordings like the baseline do show the generator's gap. The population median is low because most recordings sit elsewhere.
- **Below 0.215:** even baseline-like recordings show a small gap. The cause is in what the generator leaves out.
- **CI spans 0.215:** with n 13, step 1 is then unresolved and step 2 decides.

**Step 2, twins (`src/sim_gap_twins.py`).** Each TD recording is rebuilt from its own knee-mode fit and refit in both modes with `fit_one`:

| cell | aperiodic | peaks | noise |
|---|---|---|---|
| M0 | median (b₀, k₀, χ₀) | none | white, SD 0.02 |
| T00 | own (b, k, χ) | none | white, SD 0.02 |
| T10 | own | own knee-fit peak model | white, SD 0.02 |
| T01 | own | none | own knee-model residual (log10 PSD − aperiodic − peaks) |
| T11 | own | own peaks | own residual: the real spectrum |

White noise for recording j (features order) uses seed 20261029 + j: the same noise vectors as (a)'s replicates 0–184 at nlv 0.02. Closure check: T11 must reproduce the real +0.110.

Effects on median d. The first three add up exactly to M0 − T11:
- **Heterogeneity** of aperiodic parameters = M0 − T00.
- **Real peaks** = mean of (T00 − T10) and (T01 − T11).
- **Real residual** instead of white noise = mean of (T00 − T01) and (T10 − T11).
- **Interaction** = (T00 − T10) − (T01 − T11). Reported, not part of the sum.

Readings, as shares of M0 − T11:
- A factor **carries** the gap at 50% or more, **contributes** at 20–50%, and is **negligible** below 20%.
- If |interaction| is 20% or more of M0 − T11, peaks and noise are not separable, and only their sum is reported as identified.

Bootstrap CIs are over recordings (10,000 resamples, seed 20261029). Every cell uses the same resamples, so the effects are paired.

**Scope.**
- The twins inherit the knee-mode decomposition: its aperiodic curve, its peaks and its residual. They explain why the generator missed the real gap; they do not say which aperiodic model is right.
- Flagged recordings are included. A knee ≤ 0 still gives a valid curve on 1–49.5 Hz (k > −1 for all nine); their knee frequency is undefined, so they form their own bin.
- Per-recording twin values stay in `data/`. Aggregates go in `tables/sim_gap_twins.csv` and `tables/sim_gap_by_knee_freq.csv`.
- No change to `simulate.py` or `sim_summary.csv`.

## 2026-11-16 — Why the synthetic baseline overpredicted the TD gap: results

**Run.** `src/sim_gap_twins.py` fits 185 TD recordings × 5 cells × 2 modes, 1,850 fits, all converged, in 2 min 32 s.
- **Closure.** T11 rebuilds each recording as aperiodic + peaks + residual. It reproduces the real TD median d exactly (0.110248), with a per-recording maximum |difference| of 1.5e-5.
- **M0 against (a).** M0 matches (a)'s replicates 0–184 to within 1.8e-5 per replicate; the small difference is because b₀, k₀ and χ₀ come from `features.csv` rather than a refit. Replicate 163, the one knee failure in (a), converges here.

**Step 1** (`tables/sim_gap_by_knee_freq.csv`). Median d [95% CI]:

| knee frequency | n | real d | T00 twin d (own aperiodic only) | real fixed / knee periodic delta |
|---|---|---|---|---|
| < 2 Hz | 43 | +0.068 [0.014, 0.102] | +0.107 | 0.357 / 0.251 |
| 2–3 Hz | 25 | +0.120 [0.054, 0.198] | +0.267 | 0.205 / 0.008 |
| 3–4.2 Hz | 33 | +0.209 [0.106, 0.281] | +0.337 | 0.262 / 0.003 |
| 4.2–6 Hz | 47 | +0.158 [0.111, 0.221] | +0.383 | 0.202 / 0.006 |
| ≥ 6 Hz | 28 | +0.104 [0.062, 0.146] | +0.400 | 0.116 / 0.000 |
| no knee frequency | 9 | −0.010 [−0.147, −0.000] | −0.008 | 0.435 / 0.582 |
| near baseline (3.0–4.2 Hz, χ 2.0–2.4) | 13 | +0.106 [0.028, 0.281] | +0.337 | 0.110 / 0.000 |

The 3–4.2 Hz and near-baseline rows share their T00 median and their upper CI limit by coincidence: the same recording is the median twin in both subsets.

Under the declared reading, step 1 is **unresolved**: the near-baseline median is +0.106, but its CI spans 0.215 (n 13). The point estimate agrees with step 2 below.

Not declared, descriptive only:
- Real d peaks in the 3–4.2 Hz bin and falls on both sides.
- The aperiodic-only twin instead rises with knee frequency. Above 4.2 Hz (75 recordings) it overpredicts real d by 0.23–0.30; below 2 Hz, by 0.04.
- In the < 2 Hz bin, knee mode itself puts 0.25 of periodic delta in 2–4 Hz (median).

**Step 2** (`tables/sim_gap_twins.csv`). Median d [95% CI], with median periodic delta in each mode:

| cell | d | fixed | knee |
|---|---|---|---|
| M0 median aperiodic, white noise | +0.321 [0.314, 0.323] | 0.332 | 0.022 |
| T00 own aperiodic, white noise | +0.290 [0.267, 0.319] | 0.327 | 0.024 |
| T10 + own peaks, white noise | +0.186 [0.146, 0.212] | 0.262 | 0.005 |
| T01 own aperiodic + own residual | +0.195 [0.178, 0.209] | 0.247 | 0.031 |
| T11 the real spectrum | +0.110 [0.094, 0.129] | 0.228 | 0.015 |

Effects, as shares of M0 − T11 = +0.210 [0.189, 0.227]:

| effect | size [95% CI] | share | reading |
|---|---|---|---|
| heterogeneity, M0 − T00 | +0.031 [−0.002, +0.054] | 15% | negligible at the population median (the CI reaches 26%) |
| real peaks | +0.094 [0.078, 0.115] | 45% | contributes |
| real residual instead of white noise | +0.085 [0.062, 0.111] | 41% | contributes |
| interaction | +0.019 [−0.004, +0.055] | 9% | below 20%: separable (the CI reaches 26%) |

**Verdict under the declared readings.**
- No single factor carries the gap.
- The recordings' own peaks and their real residual contribute about equally (45% and 41%), and the two are separable. Together they account for 85% (0.18 of 0.21).
- Spread of the aperiodic parameters is negligible at the population median (15%). Its CI reaches 26%, so "contributes" is not excluded. Within bins it is not small: the T00 twin d runs from 0.11 (< 2 Hz) to 0.40 (≥ 6 Hz); the spread largely cancels at the median.
- The interaction's CI also reaches 26%, so the peaks/residual split is separable by the declared rule on the point estimate, not decisively.
- The generator overpredicted mainly because it has no peaks and only white noise at SD 0.02, not because it used one median curve.

**Where the gap goes** (descriptive, not declared; the split threshold of 0.05 was chosen after seeing the cells; `split` rows of `tables/sim_gap_twins.csv`). Median d is not the difference of the two modes' medians: from T00 to T11, the medians give 0.303 → 0.213, but median d falls 0.290 → 0.110. Paired per recording, T00 → T11 (the real spectrum) follows two routes:

| recordings | n | change in knee periodic delta | change in fixed | change in d | d: T00 → T11 |
|---|---|---|---|---|---|
| knee mode gains > 0.05 | 64 (35%) | +0.239 | +0.217 | −0.046 | 0.154 → 0.090 |
| the rest | 121 (65%) | −0.010 | −0.240 | −0.222 | 0.364 → 0.121 |

(medians). In two thirds of recordings, knee mode stays near 0 and fixed mode loses about 0.24 of 2–4 Hz periodic delta once the recording's own peaks and residual are present: this is where the gap at the median goes. In the other third, both modes pick up a real 2–4 Hz peak almost equally, so d barely moves; this is what lifts knee mode's upper tail (q75 0.04 → 0.19, q90 0.05 → 0.40) and, on the mean, outweighs the fixed-mode fall (mean change: knee +0.092, fixed −0.054). How specparam produces the fixed-mode fall (the relative peak threshold, the robust aperiodic fit) is not tested here.

**For the 2026-11-14 verdict.**
- "Size not held" stands.
- The PROTOCOL size came from a generator without peaks and with white noise; those two omissions account for most of the 3× overprediction.
- The direction clauses do not depend on the size.
- The AS synthetic values would be overstated the same way if the same omissions matter there. Not checked.

**Not done.** The residual effect mixes size and structure. The real knee-model fit error is larger than the generator's: median 0.029 against 0.011 for (a)'s knee fits at nlv 0.02. A sixth cell, with white noise at each recording's own residual SD, would separate the two. Not run.

**Reproduction.** I made a fresh copy of the tracked and new files (no `data/`), built a new venv from `requirements.txt` and ran `bash scripts/run_all.sh`, which now includes the twins step.
- **Run.** Exit 0 in 3 min 15 s, with 29/29 tests passing. This is the second such run, made after the per-recording split columns were added to `tables/sim_gap_twins.csv`; the first, before that change, also matched on all 54 outputs.
- **Outputs.** All 54 compared outputs are byte-identical to the project's: 19 tables, 4 figures and 31 local `data/` files. The arrays in `data/fits.npz` are equal.
- **Not regenerated.** As on 2026-11-14, `sim_summary.csv` and the synthetic figure were copies, since the run did not use `--with-sim`.
- **Not pinned.** The files are uncommitted, so this run is not tied to a commit.

## 2026-11-16 — Residual effect split into size and structure: design and readings, declared before the run

Post hoc and exploratory. A follow-up to the twins entry above, requested after its results. The quantity to split is the residual effect, +0.085 [0.062, 0.111]: the change in median d from replacing white noise (SD 0.02) with each recording's real knee-model residual.

**Cells.** Two new cells, not one: the residual effect is averaged over both peak levels, so the split needs a cell at each level for the parts to add up exactly.

| cell | aperiodic | peaks | noise |
|---|---|---|---|
| T0s | own | none | white, rescaled to the recording's own residual SD |
| T1s | own | own knee-fit peaks | white, rescaled to the recording's own residual SD |

The noise is T00's vector for that recording (seed 20261029 + j), rescaled so that its sample SD equals the SD of the recording's knee-model residual across the 98 bins. The existing cells are untouched. TD residual SD, 10th–90th percentile: 0.025–0.063, median 0.038. That is 1.9 times the generator's 0.02, and 96% of recordings are above 0.02.

**Effects.** These two add up exactly to the residual effect:
- **Size** = mean of (T00 − T0s) and (T10 − T1s): raising white noise from SD 0.02 to the recording's own SD.
- **Structure** = mean of (T0s − T01) and (T1s − T11): replacing white noise of the right SD with the real residual. "Structure" is everything about the real residual except its overall SD: correlation between bins, the shape of its distribution, and where along the frequency axis it is large. The last item is size within the delta band, so a "structure" result does not by itself mean autocorrelation.

**Readings,** as shares of the residual effect, with paired bootstrap CIs (same resamples as before):
- **Thresholds.** A part carries the residual effect at 50% or more, contributes at 20–50%, and is negligible below 20%. The denominator is small (0.085), so the shares will be noisy; the sizes and their CIs are the primary report.
- **Negative size effect.** If larger white noise widens the gap rather than narrowing it, structure's share exceeds 100%. Read that as: structure carries the residual effect, and size works against it. Both are reported.
- **Peak levels disagree.** If the size effect without peaks (T00 − T0s) and with peaks (T10 − T1s) have opposite signs and each is at least 20% of the residual effect in size, the split depends on peaks. It is then not identified as a single number, and both are reported.
- **Fit error.** Median knee-fit error in T0s and T1s is reported against the real 0.029. Knee mode absorbs some white noise into peaks (11 in (a)), so if the fit errors differ, the size is matched on residual SD, not on fit error.

**Already seen, before this run.** The synthetic grid varied white-noise size at the median curve: (a) gave +0.308 at nlv 0.05 against +0.321 at nlv 0.02. A change of −0.013, about 15% of 0.085, from a larger increase in SD than the median recording's. This points to structure before the run, and is recorded so the reading is not adjusted afterward.

**Outputs.** New rows in `tables/sim_gap_twins.csv`. The existing cell, effect and split rows must come back identical. Per-recording SDs and fits stay in `data/`.

## 2026-11-16 — Residual effect split into size and structure: results

**Run.** `src/sim_gap_twins.py` now has seven cells, 2,590 fits, all converged, in 4 min 16 s. The existing rows came back unchanged:
- all 12 earlier rows of `tables/sim_gap_twins.csv`, every column exactly equal;
- `tables/sim_gap_by_knee_freq.csv`, byte-identical;
- all 925 earlier per-recording d values, identical.

**Cells.** Median d [95% CI], with median periodic delta in each mode and median knee-fit error:

| cell | d | fixed | knee | knee-fit error |
|---|---|---|---|---|
| T00 own aperiodic, white 0.02 | +0.290 [0.267, 0.319] | 0.327 | 0.024 | 0.011 |
| T0s own aperiodic, white at own SD | +0.264 [0.230, 0.293] | 0.315 | 0.042 | 0.022 |
| T01 own aperiodic, own residual | +0.195 [0.178, 0.209] | 0.247 | 0.031 | 0.019 |
| T10 + own peaks, white 0.02 | +0.186 [0.146, 0.212] | 0.262 | 0.005 | 0.016 |
| T1s + own peaks, white at own SD | +0.178 [0.143, 0.198] | 0.264 | 0.024 | 0.026 |
| T11 the real spectrum | +0.110 [0.094, 0.129] | 0.228 | 0.015 | 0.029 |

**Effects,** as shares of the residual effect (+0.085) and of the whole gap (+0.210):

| effect | size [95% CI] | of residual | of whole gap | reading |
|---|---|---|---|---|
| residual size | +0.017 [+0.006, +0.040] | 20% | 8% | contributes, exactly at the threshold (20.1%) |
| residual structure | +0.068 [+0.040, +0.087] | 80% | 32% | carries |
| size, without peaks (T00 − T0s) | +0.026 [+0.012, +0.058] | 30% | | |
| size, with peaks (T10 − T1s) | +0.008 [−0.009, +0.033] | 10% | | |

**Verdict under the declared readings.**
- **Structure carries the residual effect** (80%).
- **Size contributes** by the declared rule, but only just: its point share is 20.1%, and its CI covers 6%–47% of the residual effect. Its point share sits on the threshold, so the label decides little; its CI excludes zero.
- **The peak levels agree in sign,** so the disagreement rule is not triggered. Size matters less once the recording's peaks are present (10% against 30%).
- **The prior from the grid held.** It pointed to structure, and the reading was not adjusted.
- **Fit error.** The match is on residual SD, not on fit error. With peaks, the size-matched cell's knee-fit error is 0.026 against the real 0.029; it covers most of the step up from white 0.02 (0.016). Knee mode absorbs some white noise into peaks, which accounts for the rest.
- **What "structure" covers.** As declared, it includes where along the frequency axis the residual is large, so this result does not by itself mean autocorrelation. It was not split further.

**The whole gap, combined with the twins entry.**

| what the generator left out | share of the 0.21 |
|---|---|
| the recording's own peaks | 45% |
| the structure of the real residual | 32% |
| the recording's own curve, instead of the median one | 15% |
| the size of the real residual | 8% |

The four add up to 100%. The peaks × residual interaction (9%) is reported separately. The generator's 3× overprediction comes mainly from the absence of real peaks and from white noise standing in for a structured residual; the noise level and the median-curve shortcut matter little at the median.

Not done: a generator that resampled real residuals and real peaks would test whether its predicted sizes then match. This is a suggestion only.

**Reproduction.** I made a fresh copy of the tracked and new files (no `data/`), built a new venv from `requirements.txt` and ran `bash scripts/run_all.sh` with the seven-cell twins.
- **Run.** Exit 0 in 5 min 0 s, with 31/31 tests passing.
- **Outputs.** All 54 compared outputs are byte-identical to the project's: 19 tables, 4 figures and 31 local `data/` files. The arrays in `data/fits.npz` are equal.
- **Not regenerated.** `sim_summary.csv` and the synthetic figure were copies, since the run did not use `--with-sim`.
- **Not pinned.** The files are uncommitted, so this run is not tied to a commit.

## 2026-11-16 — Phase 4 complete

**What Phase 4 covered.**
- *2026-11-14:* the primary analysis, both predeclared sensitivities, secondaries (a)–(e) and the AS–TD contrast by mode.
- *2026-11-14, post hoc:* the periodic-gap subsets and Fig 1e cross-tabs.
- *2026-11-16, post hoc:* why the synthetic validation overpredicted the TD gap.

**Headline.**
- **Primary median d** (fixed − knee periodic delta): all +0.197 [0.161, 0.229], 92% > 0; TD +0.110 [0.094, 0.129]; AS +0.426 [0.349, 0.481].
- **Sensitivities:** both agree.
- **AS–TD contrast:** it grows under fixed mode. Cohen's d for periodic delta goes from 1.27 to 1.72, and for exponent from 0.87 to 1.53.
- **TD exponent age slope:** +0.320 under knee against +0.078 under fixed, per doubling of age.

**Predeclared item not done.** The GAM shape comparison in secondary (c), against the paper's Fig 3 TD curve, was not run, and no reason was recorded. The "can say" question about the TD curve's shape is therefore answered for the linear slope only.

**Records written at the close.**
- **`README.md`** has a Results section: the primary with its CI, both sensitivities, the contrast under each mode, Fig 1e in two sentences, and PROTOCOL v2.3's "It cannot say" list copied verbatim, including the upstream-processing (aggregation-level) item.
  - That list replaces the README's earlier "What this study cannot claim" section. Its third bullet said the reported values "are averages of regional fits", which is more than the Fig 1e diagnostic found: Fig 1e is not a single knee fit, and which aggregation was used is unknown.
  - The Status line, run line and Layout were updated to match the pipeline.
- **`tables/deviations.csv`** has 11 new rows of differences from PROTOCOL, each with its reason: 3 procedure, 1 specification, 1 amendment, 3 expectation, 2 addition and 1 exploratory plan.
  - The 7 reported-value rows are unchanged on their original columns.
  - The rows live in `src/validation.py` (`PROTOCOL_DEVIATIONS`), so `run_all.sh` regenerates them. Each reason is what NOTES recorded at the time, or the row says none was recorded.
- **`PROTOCOL.md`** is unchanged (v2.3).

**Reproduction.** I made a fresh copy of the tracked and new files (no `data/`), built a new venv from `requirements.txt` and ran `bash scripts/run_all.sh`.
- **Run.** Exit 0 in 5 min 0 s, with 31/31 tests passing.
- **Outputs.** All 54 compared outputs are byte-identical to the project's: 19 tables including `deviations.csv`, 4 figures and 31 local `data/` files. The arrays in `data/fits.npz` are equal.

**State.**
- **Uncommitted.** Every Phase 4 file is uncommitted on top of e0f3d0a, so this reproduction is not tied to a commit.
- **Synthetic outputs.** `sim_summary.csv` and `fig_synthetic_validation.png` were last regenerated and verified at 7a2bf66 (2026-11-14 gate entry).

**Still open.**
- The AS periodic-delta gap against the reported values (−0.125) is unexplained. How the authors aggregated is the email's question.
- The mechanism clause of the 2026-11-03 expectation is untested.
- The GAM shape comparison has not been run.
- No synthetic generator that resamples real peaks and residuals has been built.

## 2026-11-19 — PROTOCOL v2.4: post-analysis records

`PROTOCOL.md` moves from v2.3 to v2.4. It records the Phase 4 close by adding to the document; nothing is rewritten.
- **Diff.** The only line removed is the title. Every addition is marked as written 2026-11-19 and cites the dates of the results it reports.
- **Predeclared text.** Unchanged word for word: the estimand, secondary outcomes (a)–(e) including the GAM clause in (c), inclusion, conventions, the QC gate, statistics ("log age"), the expectations and the claims section.

**Added.**
- **Under the title:** a note that predeclared content froze at v2.3 (2026-11-05) and that v2.4 adds post-analysis records only.
- **After the Conventions deviations sentence:** a pointer to the protocol rows of `tables/deviations.csv`.
- **Expectations:**
  - a closing line to the real-data exposure note: nothing predeclared remains unseen as of 2026-11-14;
  - the outcome of each clause, with numbers from `tables/`. The TD size was not held. The mechanism is untested and therefore not confirmed. The reason for the TD miss is labelled post hoc.
- **Exploratory:** the four listed analyses were not run, and the exploratory and post hoc analyses that were run are listed instead. The AS − TD median difference appears only here; it is not added to the secondary outcomes.
- **New section, Deviations from this protocol:** placed after the claims section and before the change log. It points to `tables/deviations.csv` and lists the items that touch predeclared content:
  - the GAM in (c) was not run, and no reason was recorded;
  - the committed paired plot is binned;
  - (b) is reported in `group_contrast_by_mode.csv`;
  - "log age" was taken as log2;
  - QC item 2 was resolved by the v2.3 rewording.
- **Change log:** a v2.4 entry saying it was written after all results were seen.

**Elsewhere.**
- **README:** the copied "cannot say" list now cites "v2.3; unchanged in v2.4", and it still matches PROTOCOL word for word.
- **`src/validation.py` and `tables/deviations.csv`:** unchanged. They name the differences "from PROTOCOL v2.3", which still describes the predeclared plan correctly. No output was regenerated, because no file read by the pipeline changed.

Nothing committed.
