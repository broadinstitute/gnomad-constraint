# Ridge-penalized rr: summary

Condensed from two documents.

- **[Follow-up report][doc]** (`rr_syn_dnm_followup.md`, dated Nov 5, branch `jg/rr_ridge_fit`): the analysis behind the choice.
- **[Slides, "Constraint code reorg and the regional mutation rate"][art]** (Nov 17): the implementation plan.

**rr** is the regional factor that scales expected counts; our deck calls it adj_r. Gnocchi's released table is `adj_r`, and the legacy synonymous refit is `syn_adj_r`.

## Decision

- **Model:** fit rr with a **ridge penalty**. This replaces both the Bonferroni feature screen and the strength α.
- **Training data:** the **Satterstrom 2026** synonymous DNMs, against the current background. No alternative background did clearly better.
- **Normalize per context,** so that rr averages 1 over exome synonymous sites.
  - Without this, Gnocchi `adj_r` put synonymous o/e about 8–10% low.
  - Raw synonymous rr adds about 1.4% to every expected count, LoF included.
- **No rr on chrX/Y.** Six of the 13 features have no X/Y data, and the best X/Y option lowers noise by under 1%.

## How the ridge rr works ([details][model])

- **One model per context.** Each of the 32 trinucleotide contexts gets its own logistic regression:
  - outcome 1: synonymous sites with a DNM;
  - outcome 0: a 5% sample of possible synonymous sites with no DNM and no variant in gnomAD v4 exomes, one row per locus.
- **Features.** The 13 Gnocchi regional features are each measured at 4 scales (1 kb, 10 kb, 100 kb and 1 Mb), giving 52 columns, z-scored within the context. CpG contexts drop the 5 features linked to methylation.
- **The penalty.** Ridge keeps every column and shrinks every weight toward 0. It minimizes mean log loss + (λ/2)·Σβ², with the intercept left unpenalized.
- **Choosing λ.** λ = 0.3 per training row, the value with the lowest DNM log loss on held-out chromosome halves (searched over 0.05–1).
  - The search uses no gnomAD data, so the checks on gnomAD stay independent.
  - Retune λ if the background or the model changes ([penalty][pen]).
- **Turning the model into rr.** For each 1-kb window and context, rr = σ(β₀ + β·z) / σ(β₀) ≈ exp(β·z). rr is then normalized and multiplied onto each variant's expected count, after the plateau and coverage models.
- **Size of the effects.** Each feature changes the rate by about 2% per SD or less. rr spans about 0.9–1.1 in A/T contexts and 0.85–1.25 in C and CpG contexts.

## Why ridge ([details][strength])

Without a penalty, every rr version was too strong; the Bonferroni port had a calibration slope of 0.61. There are two causes:

- features picked by p-value have inflated effects;
- the model is fit on DNMs but applied to gnomAD.

α only corrected the symptom; ridge removes the cause.

| rr version | Noise | Noise, held-out chr | Slope (1 = right strength) | Syn outliers, no QC flag | LoF expected |
|---|---|---|---|---|---|
| No correction | 2.91 | 2.91 | — | 123 | — |
| Bonferroni screen, α 1 (the port) | 2.65 | 2.70 | 0.61 | 79 | −2.5% |
| 1e-2 screen + α 0.625 (October pick) | 2.49 | 2.50 | 0.82 | 74 | not run |
| **Ridge** | **2.45** | **2.47** | **0.86** | **71** | **−1.1%** |
| Roulette + ridge (one untuned try) | 2.06 | 2.09 | 0.77 | 37 | −5.1% |

How to read the table:

- **Noise** is the per-gene synonymous dispersion: mean z² on the core gene set. Chance alone gives at most 0.70.
- **The ridge row** is the run with QC-flagged genes left out of training. With all genes, ridge's noise is 2.46, or 2.48 on held-out chromosomes.
- **Gnocchi `adj_r`** (not in the table) is worse than no correction: noise 5.45, slope 0.26.

Other findings:

- **No α needed:** ridge's best α is 1.
- **The DNM set stops mattering.** The gap between the old and new sets falls from 0.155 to 0.004. The new set only looked worse because a fixed p-value cutoff lets in more weak features as the DNM count grows ([details][dnm]).
- **Disease genes are unchanged or slightly better:**
  - LOEUF rank correlation with the no-rr run: 0.992.
  - ClinGen HI AUC: 0.867 → 0.870.
  - NDD AUPRC: 0.144 → 0.147, and 0.081 → 0.088 for genes with under 50 expected pLoF.
  - HI AUPRC: 0.168 → 0.165, within its interval.
- **Missense:** expected counts move by −0.4%, and 9% more genes reach missense z ≥ 3.09.

## Implementation plan (slides)

The port (`rp/add_mu_adjustment`: +4,398 lines in 12 files, based on an April commit) is being rewritten as small PRs. They sit on top of a code reorg that only moves code.

| PR | What it does | Status (GitHub, Nov 17) |
|---|---|---|
| gnomad_methods #862 | Splits `constraint.py` into the `gnomad/utils/constraint/` package; old imports still work | merged Nov 13 |
| gnomad_methods #863 | Ridge fit in Hail (Newton's method): `regional_mutation_rate.py` and tests | open |
| gnomad-constraint #51 | Per-base pipeline; the bottom of the stack | open |
| #67 | Splits `constraint_pipeline.py` into `prepare_data.py`, `expected_models.py`, `expected_counts.py` and `constraint_metrics.py` (moves only) | open |
| #64 | Imports the DNMs and the Gnocchi feature matrix (`prepare_data.py`) | open |
| #65 | `--apply-regional-rr`: in the aggregate step of `constraint_metrics.py`, expected × rr, or × 1 where there is no rr | open |
| #66 | Training steps in `expected_models.py` (site features, training set, penalty tuning, per-context models), and `--compute-rr` in `expected_counts.py` | open |
| #57 | Generated flowchart page; a repo admin has to turn on GitHub Pages | open |

The gnomad-constraint PRs stack in the order #51 → #67 → #64 → #65 → #66 → #57.

What the PRs do:

- **Settings:** λ = 0.3 per site, from #66's tuning step on held-out DNMs; training on the Satterstrom 2026 ASD and DDD files; rr = 1 on chrX/Y and in windows without features.
- **When rr is applied:** after calibration. The plateau and coverage models aren't refit.
- **Normalization:** `--compute-rr` normalizes rr on the per-variant expected Table. `--apply-regional-rr` refuses an rr Table that was normalized on different models.
- **Release:** records the settings in `regional_rr_params`.
- **Opt-in:** rr is off unless `--apply-regional-rr` is passed. The port applied it by default.

## To check before merging

- **Score the PRs' own output.** Every number above comes from the `jg/rr_ridge_fit` runs. The #863 Hail fit was checked only on the old training Table: all 576 fits converge, and tuning picks 0.3.
- **Training DNMs.** The PRs import only the Satterstrom files, which keeps 20,552 of the 23,336 DNMs the full run trained on. Confirm which set is intended.
- **QC-flagged genes in training.** The report recommends leaving them out (a small gain), and its ridge numbers come from that run. The slides don't mention it.
- **Run end to end with rr.** `--compute-constraint-metrics` and `--prepare-release` haven't been run with rr yet. They need `hl.qgamma` (Hail ≥ 0.2.137), which the test cluster's Hail lacks.
- **The gnomad_methods pin.** #66 pins #863's head (`8bb46c58`); move the pin to #863's merge commit before #66 merges.
- **Items the report raised on the port.** The slides don't say whether the rewrite settles them:
  - the `unadj_*` release fields;
  - the corrected 1-bp window coordinates;
  - #55's sums of the non-coding `adj_r`;
  - the NA-observed fix (#2298), which needs the context Table rebuilt.

## Open decisions (slides)

- **Roulette's per-variant rates: worth a proper look?** One untuned try, with ridge on top ([details][rates]):
  - **For:** it removed about twice as much noise as ridge alone (2.06 against 2.45), fixed the 8p maternal hotspot (o/e 1.02, against 1.21 with ridge), and also won on held-out DNMs (+1,316 log-likelihood units over chance, against +114 for ridge).
  - **Against:** the combination is too strong (slope 0.77), and it moves LoF by −5% and missense by +3% relative to synonymous, which the calibration can't check. Part of its gain on gnomAD is gnomAD fitting gnomAD.
- **chrX/Y:** keep rr at 1, or rebuild the missing features?
- **Default:** turn `--apply-regional-rr` on by default for the next release?
- **gnomad_methods:** release a version after #863 instead of pinning commits?

## What ridge leaves unfixed (report)

- **Leftover noise** ([details][noise]):
  - About half comes from individual sites, which is a job for the mutation model.
  - The other half is gene-wide. DNMs confirm about 70% of that gene-wide spread as real rate differences.
- **Strength differs by context class.**
  - CpG: far too strong in every version (slope 0.10–0.25).
  - A/T: too weak with ridge (slope about 2), and the DNMs show no trend with rr at all.
- **Outliers** ([details][outliers]):
  - About 130 genes with far fewer variants than expected are mapping artifacts (`low_exome_mapping_quality`); leave them to the flags.
  - The chr6p22 histone cluster (24 outliers, o/e about 1.85) is gene-specific, not regional, so no rr can fix it. Flag it, or adjust it gene by gene.
  - Half of the other outliers with excess variants (24 of 46) sit in maternal mutation hotspots, which ridge only partly corrects (8p: 1.32 → 1.21).
- **Overall synonymous o/e of 1.015.** This comes from the plateau calibration, not from rr. A saturating plateau curve brings it to 1.000.
- **Later ideas** ([details][later]):
  - one model for expected counts in place of the stack (noise 1.65, against 2.43 for the stack with ridge);
  - one rr strength per context class;
  - a gene-level adjustment, shrunk toward 1.

[doc]: https://github.com/broadinstitute/gnomad-constraint/blob/jg/rr_ridge_fit/gnomad_constraint/plots/adj_r_syn_oe/rr_syn_dnm_followup.md
[art]: https://claude.ai/code/artifact/875ba016-0900-4f4f-8fde-70a2eb377209
[model]: https://github.com/broadinstitute/gnomad-constraint/blob/jg/rr_ridge_fit/gnomad_constraint/plots/adj_r_syn_oe/rr_syn_dnm_followup.md#how-the-rr-model-works
[pen]: https://github.com/broadinstitute/gnomad-constraint/blob/jg/rr_ridge_fit/gnomad_constraint/plots/adj_r_syn_oe/rr_syn_dnm_followup.md#how-each-version-uses-the-features
[strength]: https://github.com/broadinstitute/gnomad-constraint/blob/jg/rr_ridge_fit/gnomad_constraint/plots/adj_r_syn_oe/rr_syn_dnm_followup.md#strength-details
[dnm]: https://github.com/broadinstitute/gnomad-constraint/blob/jg/rr_ridge_fit/gnomad_constraint/plots/adj_r_syn_oe/rr_syn_dnm_followup.md#dnm-sets-details
[rates]: https://github.com/broadinstitute/gnomad-constraint/blob/jg/rr_ridge_fit/gnomad_constraint/plots/adj_r_syn_oe/rr_syn_dnm_followup.md#per-variant-rates-roulette-and-how-the-pipeline-uses-it
[noise]: https://github.com/broadinstitute/gnomad-constraint/blob/jg/rr_ridge_fit/gnomad_constraint/plots/adj_r_syn_oe/rr_syn_dnm_followup.md#what-the-leftover-noise-is
[outliers]: https://github.com/broadinstitute/gnomad-constraint/blob/jg/rr_ridge_fit/gnomad_constraint/plots/adj_r_syn_oe/rr_syn_dnm_followup.md#the-synonymous-outliers
[later]: https://github.com/broadinstitute/gnomad-constraint/blob/jg/rr_ridge_fit/gnomad_constraint/plots/adj_r_syn_oe/rr_syn_dnm_followup.md#still-worth-trying-or-looking-at
