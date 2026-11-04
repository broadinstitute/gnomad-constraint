# Synonymous-DNM rr: where things stand (2026-10-27)

## The problem, in plain terms

- Constraint compares how many variants we **observe** in a gene with how many we **expect** without selection.
- The expected count comes from a mutation model that only looks at each site's immediate sequence: the base, its two neighbors, and CpG methylation.
- Mutation rates also vary by **region**, with GC content, recombination and similar features. The model ignores this. Genes in more mutable regions therefore look less constrained than they are, and genes in less mutable regions look more constrained.
- Synonymous variants are the test, because they are nearly neutral and their obs/exp should be about 1 in every gene. It isn't:
  - genes scatter far more than chance allows (mean z² ≈ 2.9, where chance alone gives at most 0.7)
  - about 270 genes are extreme outliers, roughly half of them known mapping artifacts
- **The fix under test:** multiply each site's expected count by a **regional factor (rr)**, learned from where de novo mutations land. The open question is which rr to use, and how to apply it.

## Recommendation

Use a ridge-penalized rr, normalized per context, with no rr on chrX/Y. Then confirm it with one full pipeline run. It's the black "Ridge, Satterstrom" series in every figure.

1. **Normalize rr per context.** Expected counts are already calibrated to exome synonymous sites. Raw rr adds about 1.4% to every gene, LoF included, and gains nothing in return.
2. **No rr on chrX/Y.** The best option there lowers noise by under 1%, and correcting only some contexts is hard to explain. Revisit this if the missing features can be recovered.
3. **Replace the Bonferroni screen and α with a ridge-penalized fit.**
   - **How it works:** every feature goes in, and the penalty decides how much each one counts. The penalty is chosen on held-out DNMs only, with no gnomAD data.
   - **Noise:** it has the lowest noise of any rr fit to de novo mutations: 2.46, vs 2.51 for legacy and 2.56–2.70 for Ruchit's tables.
   - **No α needed:** its calibration slope is 0.86, and adding α would gain under 1%.
   - **The DNM set stops mattering:** with ridge, the two sets perform the same.
   - **Small code change:** the per-context fit currently has no penalty (`alpha=0`). The ridge penalty is per site (mean log loss + penalty / 2 × the sum of squared coefficients), so it carries over unchanged from the held-out-chromosome fits to the full fit. A penalty of 0.3 gave the lowest held-out DNM log loss for both DNM sets. This is implemented on the local branch `jg/rr_ridge_fit` (not pushed).
4. **Choose the DNM set on other grounds.** I'd lean toward the Satterstrom set: more DNMs, more on chrX, and newer annotation. Use the union with the 17k set if the 282 dropped DNMs turn out to be re-annotation losses.
5. **Confirm it in Ruchit's pipeline, then do one full run.** These results come from a local re-implementation on his data.
   - Rerun his evaluation with ridge.
   - Do one full constraint run, combined with the NA-observed fix rebuild.
   - Check LOEUF on known haploinsufficient genes.
   - Count outliers on genes without QC flags.
6. **Research, not a blocker:** fit regional effects on gnomAD rare variants (intronic, from genomes). Locally, the version fit on exome synonymous counts beat everything by about 4%.

**Not worth effort:**
- **The top-LOEUF-decile rise:** it mostly comes from how LOEUF sorts genes.
- **The ~130 deficit outliers:** they are flagged mapping artifacts.
- **Rebuilding the background:** the test changed strength but not noise.

## What a good outcome looks like

Nobody has written these down as a formal bar; they are the checks the reports and the Slack thread have used.

| # | The rr should… | Measured by | Reference: 17k DNMs, Bonferroni, α 1 | Recommended: ridge, Satterstrom DNMs |
|---|---|---|---|---|
| 1 | Not shift the overall level that the model calibration sets | mean rr over exome synonymous sites = 1 | Yes if normalized (mean rr is exactly 1). Raw rr raises every gene's expected count by about 1.4% | Yes (normalized) |
| 2 | Reduce per-gene noise | mean z² vs. no correction (panel a) | Yes: 2.91 → 2.56 | Yes: 2.91 → 2.46 |
| 3 | Remove synonymous outliers without creating new ones | outlier counts (panel b) | Yes: 270 → 218 (61 removed, 9 new) | Yes: 270 → 218 (unflagged 123 → 71) |
| 4 | Not overcorrect: genes rr calls more mutable really do carry more variants | calibration slope, 1 = perfect (panel c) | Acceptable: 0.68. Gnocchi `adj_r`, at 0.28, fails | Better: 0.86; adding α would gain under 1% |
| 5 | Hold up on chromosomes it wasn't trained on | held-out-chromosome noise | Yes: 2.61 vs. 2.91 | Yes: 2.48 vs. 2.91 |
| 6 | Leave LoF and missense metrics sensible: LOEUF still separates known haploinsufficient genes, flags look reasonable | full pipeline run, compared with the unadjusted values | **Not done yet** | **Not done yet** |
| 7 | Be simple to explain and maintain: pipeline code, close to the legacy method, one consistent rule for chrX/Y | — | Mostly; **the chrX/Y rule is still open** | Yes: no p-value cutoff, no α, no rr on chrX/Y |

Sources:
- Ruchit's `rr_syn_dnm_report.html` (v1, 10/16) and `rr_syn_dnm_report.v2.html` (10/19–21)
- the Slack thread
- your 10/1 report, `investigation/adj_r_investigation_report.md` on `jg/adj_r_syn_oe_analysis`
- `origin/rp/add_mu_adjustment` at `53a8f1d`

## What you need to do

1. **Answer Ruchit's open question ("normalize rr per context?"): yes.**
   - Normalization was the fix your report proposed. His data show it changes only the overall level, not the noise. Level-free dispersion for legacy is 4.378 raw vs 4.373 normalized.
   - Raw rr averages about 1.014 over exome synonymous sites, so it raises every gene's expected count by about 1.4%, LoF and missense included. That shifts LOEUF for every gene.
   - Don't argue this from syn o/e by decile. Whether the constrained deciles sit at 1.00 or 1.01 depends on which flagged genes are dropped (panel g), and normalization only moves the whole curve by about 1.3%.
   - Raw rr only looks better because uncorrected syn o/e is 1.015 overall, driven by the least constrained genes. That offset belongs in the model calibration, not in rr.

2. **Answer Kaitlin ("legacy didn't have a strength adjustment").** She's right. The strength α is new, and it came from your report.
   - Your report called legacy `syn_adj_r` "well calibrated" (slope 0.70, best α 0.75). Ruchit calls the same numbers "25–30% too strong."
   - At the Bonferroni cutoff, fitting α gains little: legacy core-set dispersion goes from 2.509 to 2.461 (about 2%), and it costs 6–8 synonymous outliers. α matters only with the Satterstrom DNMs or a looser cutoff.
   - With the old DNM set, dropping α (Ruchit's latest pick) is consistent with your report. With the Satterstrom set, α matters more (best 0.625, 7% less noise). So the strength decision waits on the training set.
   - kchao's request for a per-variant version of the strength figure is already met: v2 section 9 applies α per site and gets the same α (0.75 for legacy).

3. **Tell Ruchit that the ~130 "deficit" outliers are flagged mapping-artifact genes** (new check, figure at the end).
   - v2 says about 130 deficit synonymous outliers (obs ≪ exp) are unchanged by every rr option, and lists characterizing them as a next step.
   - On his own evaluation set, 129 of the 130 carry the release flag `low_exome_mapping_quality`, with median syn o/e about 0.2. These are the segdup and paralog genes from your 9/28 session (RIMBP3, NBPF1, …).
   - His outlier counts use all 18,406 transcripts, although his dispersion core set already drops flagged genes. Counting unflagged genes only, syn outliers fall from 123 to 66 with legacy `syn_adj_r` and to 71–76 with Ruchit's tables, about half. All transcripts gives 270 → 217–219.
   - Suggest he count outliers on unflagged genes. That also makes the options easier to tell apart.

4. **Pin down the chrX/Y decision.** Ruchit's "no masking" is not the same as kchao's "no correction."
   - Unmasked means a partial correction. On X/Y, a context gets rr only if none of its selected features is missing there. With the old DNM set that is 3 of the 13 modeled contexts (ACT, CAC, TCT), which is the "handful of contexts" kchao objected to.
   - The benefit is tiny. Over the 769 X/Y transcripts, level-free dispersion goes from 4.364 to 4.332 unmasked, and masking makes it worse (4.435). So no rr on X/Y costs almost nothing and is the easiest to explain.
   - On kchao's "why only 6/13": the six features with no X/Y data are `recomb_male`, `recomb_female`, `met_sperm`, `Nucleosome`, `cDNM_maternal_05M` and `cDNM_paternal_05M`. They come from `RR_XY_MASKED_FEATURES` on his branch.
   - Your earlier finding explains why the gap exists: Gnocchi's published 13-feature matrix covers chr1–22 only. Ruchit rebuilt the X/Y features from the source tracks (`build_xy_features.py`), and these six have no X/Y source data.

5. **Get ready to review `rp/add_mu_adjustment`. It touches your code (#51, #55), and there's no PR yet.**
   - It is based on `5bc7d3a`, 22 commits behind `jg/switch_constraint_to_per_base` (#51), so it needs a rebase onto #51, including the new tests. The release gains `unadj_*` fields, so `tests/data/constraint_metrics_row_type.txt` and `tests/data/release_tsv_columns.txt` have to change.
   - rr is **on by default** in `--aggregate-per-variant-expected` (`--skip-regional-rr` turns it off). The release keeps the unadjusted `exp`, `oe`, `oe_ci`, `z_raw` and `z_score` alongside the adjusted values. That is a release-schema decision.
   - It overlaps #55. Both add `additional_fields_to_sum`, and #55 also sums the non-coding `adj_r`, which every analysis now says not to use for coding. Decide whether those sums come out of #55.
   - The branch now trains on the Satterstrom 2026 DNMs, which replaced the old DDD+ASD file (`644eca7`). Until the training set is chosen, both should be selectable rather than one replacing the other.
   - Your 1 bp window offset now has a `--corrected-window-coords` flag. It is off by default, which reproduces legacy. Decide whether the production table uses the corrected windows.
   - Plan a single full rerun. The NA-observed fix (#2298) needs the preprocessed context HT rebuilt, so do it in the same run as the final rr.

6. **One idea worth raising (my hypothesis).**
   - The rr background is possible synonymous sites that are absent from gnomAD v4.1 exomes (`create_rr_training_ht`, the same design as the legacy notebook).
   - High-mutation windows have more polymorphic sites, so more of their sites are removed from the background. That inflates the DNM-to-background ratio exactly where rr is high. It would make every rr too strong, which fits α ≈ 0.75, and it ties training to gnomAD, the same data rr is scored on.
   - An approximate local test supports part of it: including gnomAD-observed sites moves the slope toward 1 (0.82 → 0.93), but it doesn't lower the noise (see "Trying the alternatives").

**Not yours; just keep track of:**
- the 282 DNMs that drop out of the Satterstrom files (Ruchit and Chiara)
- Ruchit's three diagnostic tests on the larger set
- his list of the missing X/Y features

## Questions from the thread, answered

**Ruchit** = his answer. **Added** = new suggestion, not from Ruchit.

- **Why are 6 of 13 features missing on chrX/Y?**
  - **Ruchit:** their source tracks (recombination, sperm methylation, nucleosome, clustered DNMs) have no chrX/Y data. He is listing why for each.
  - **Added:** male recombination couldn't be filled in anyway, because the male X doesn't recombine outside the PARs.
- **Does correcting only some contexts on chrX/Y make sense?**
  - **Ruchit:** without masking, only 3 of 13 contexts get an rr. Masking corrects all of them but adds noise. He prefers no masking.
  - **Added:** even the best option lowers per-gene synonymous noise on chrX/Y genes by under 1% (mean z² 4.36 → 4.33, unmasked). So skip rr on chrX/Y unless the features can be recovered.
- **Normalize rr per context, or keep it raw?**
  - **Ruchit:** normalizing changes only the overall level. Raw rr would raise LoF expected counts by about 1.4% and shift every LOEUF.
  - **Added:** normalize. Expected counts are already calibrated so that exome synonymous obs/exp averages 1. rr should only move expected counts between regions, not add to the total. Raw rr adds about 1.4% on top, the same mistake (at a smaller scale) that left Gnocchi `adj_r` about 9% low. It gains nothing in return: per-gene noise is the same with or without normalization.
- **Is a strength adjustment needed? Legacy didn't have one.**
  - **Ruchit:** every rr fits best below full strength. For the 17k set the gain is only 2–3%; for the Satterstrom set it's about 7%.
  - **Added:** it depends on the DNM set because each set is too strong by a different amount. With the 17k set, full strength is only 2–3% worse than the best α and leaves about 8 fewer outliers, so α isn't worth adding. With the Satterstrom set (or a looser cutoff), full strength is clearly worse, so it needs α or a penalized fit. In a local test on held-out chromosomes, a ridge-penalized fit made α unnecessary with either set.
- **Why would more DNMs do worse?**
  - **Ruchit:** the fixed cutoff lets in weak features that predict DNMs but carry over poorly to gnomAD.
  - **Added:** features chosen by p-value have inflated effects. A penalized fit (there is none now) corrects for that. In a local test with a ridge penalty, the two sets performed the same.
- **Why do 282 old DNMs drop out, and are they still synonymous?**
  - **Ruchit:** the new files replaced the old one rather than adding to it. Most of the 282 are DD variants, which the providers re-annotated. He is confirming with them.
  - **Added:** if re-annotation is the cause, train on the union.

## Pipelines at a glance

![Shared constraint pipeline and the three ways the rr table is built](followup/rr_pipelines_diagram.png)

- **The same in every version:** steps 1–4 and 6, the per-context model, and rr defined as a window's rate divided by the average window's.
- **From adj_r to syn_adj_r:** only the training data changed, to synonymous DNMs against synonymous sites. With fewer DNMs, only 12 of the 32 contexts got a model.
- **From syn_adj_r to Ruchit's pipeline:** the same design, now built into the pipeline. It adds:
  - deduplicated DNMs, including chrX
  - chrX/Y features
  - a tunable feature-selection cutoff
  - per-context normalization
  - rr applied by default
- **Ridge fit (fourth lane; the recommended option):** Ruchit's pipeline with a ridge-penalized fit in place of the Bonferroni screen, normalized, and with no rr on chrX/Y.

### How each version performs

Black is the ridge fit with the Satterstrom DNMs (the recommended option). Dotted lines show what theory or chance predicts.

![Performance of each rr version](followup/rr_performance.png)

- **a, c:** Gnocchi `adj_r` fails on both. It nearly doubles the noise and has a slope of 0.28.
- **a:** chance alone would give a mean z² of at most 0.70, because each site is either observed or not. The best rr removes only about 18% of the noise above that level (2.91 → 2.51).
- **a, b:** legacy `syn_adj_r` and Ruchit's pipeline with the 17k DNMs are nearly tied. The Satterstrom DNMs do worse at full strength.
- **b:** by chance there would be essentially no outliers.
  - Gray is genes that already carry a release flag: about 150 under every synonymous rr, mostly the mapping-artifact genes (see the outlier figure).
  - Color is unflagged genes, the part rr changes: 123 → 66–76.
- **c, d and e measure the same thing two ways: how strong rr should be.**
  - **c** looks gene by gene. When rr says a gene should get X% more synonymous variants, gnomAD shows about slope × X% more.
  - **d** applies rr^α and finds the α with the lowest noise. α = 0 means no correction and α = 1 means rr as built. If rr were exactly right, the curve would bottom out at α = 1 (the dotted theory curve). Every rr bottoms out below 1, and that is what "too strong" means.
  - **e** shows the link: the best α in d lands on the slope from c, as theory says it should. If the real effect is slope × rr, then rr^α leaves the least error when α = slope.
  - c gives where the bottom is, and d gives how much it costs to miss it. For legacy and the 17k set, using α = 1 costs only 2–3%. Gnocchi `adj_r` bottoms out at α ≈ 0.25 and reaches 5.45 at α = 1, off the chart.
- **f:** on chrX/Y the versions differ by less than 2%, so skipping rr there costs almost nothing. The dotted line marks what rr might do on X/Y if it helped as much as on autosomes (−12%).

Panels a and c–f come from Ruchit's v2 report, sections 1, 4, 5 and 9 (v4.1.1 expected counts, 18,406 Ensembl MANE Select transcripts). Panel b is recomputed from his evaluation export. The recommended rr (Satterstrom DNMs, fit on all chromosomes) is scored the same way, with code that reproduces his numbers. Every rr is normalized except Gnocchi `adj_r`. The figures are drawn by `followup/make_followup_figures.py`.

![Synonymous obs/exp by LOEUF and s_het decile](followup/rr_by_decile.png)

- **g:** by LOEUF decile, Gnocchi `adj_r` sits about 8% low everywhere. Normalized and non-normalized `syn_adj_r` run about 1.3% apart in every decile. Every version still rises in the top deciles: to 1.11 with no correction and 1.06–1.08 with a synonymous rr.
- **h:** by s_het decile, every synonymous rr stays nearly flat at 1.00–1.05, with only a small uptick in decile 9.
- **i:** with no correction, the top decile rises to 1.11 when genes are binned by LOEUF, but only to 1.04 by s_het. Most of the rise comes from how LOEUF sorts genes, not from mutation rate.
- **j:** Ruchit's rr and the recommended rr behave the same. They remove about a third of the LOEUF rise (decile 9: 1.11 → 1.07–1.08) and stay within 1.01–1.05 by s_het.

Panels g–j are remapped from Ruchit's evaluation export (`gs://gnomad-tmp/rp/gnomad_v4.1.1/constraint/regional_rate/rr_eval_export/`), with his actual rr tables applied:
- **Data:** v4.1.1 expected counts; the same 16,241 autosomal Ensembl MANE Select transcripts in every panel, each with a v4.1.2 LOEUF decile and an s_het.
- **Filter:** genes with a release QC flag (mapping quality or coverage) are dropped, but outlier-flagged genes are kept. The outlier flags come from the same uncorrected z-scores being plotted, so dropping them would flatten the curves by construction. With his full set (QC-flagged genes kept too), the curves match his v2 section 3 exactly.
- **s_het:** gnomAD v4 s_het from Kar et al. 2026, Supplementary Table 3.
- **Checks:** the export reproduces his no-correction noise exactly (2.908, and 5.081 over all transcripts), and his normalized tables average exactly 1.
- **Code:** `followup/remap_rr_to_shet.py`, run locally.

## Remaining decisions

- **Normalize rr per context?** Recommend yes (item 1).
- **chrX/Y:** no rr, partial rr (unmasked), or recover the missing features first (item 4).
- **Training set: agreed, the Satterstrom 2026 DNMs.** Ruchit's 10/27 follow-up and this report both choose them, and the branch already defaults to them.
  - Still open: whether to take the union with the 17k set, which depends on why 282 old DNMs drop out.
- **Model fit: not settled.**
  - Ruchit recommends a 1e-2 cutoff with strength α 0.625.
  - This report recommends a ridge-penalized fit with no cutoff and no α. It scored better in a head-to-head on his data (see "Ruchit's 10/27 follow-up").
  - The ridge fit is implemented on the local branch `jg/rr_ridge_fit`.
- **For the PR:**
  - whether rr is on by default
  - the `unadj_*` release fields
  - whether the corrected 1 bp windows are used
  - whether the non-coding `adj_r` sums leave #55 (item 5)

## Still worth trying or looking at

- **Run the full pipeline once with the final rr.** Compare LOEUF, pLI and flags with the unadjusted values, and check that LOEUF still separates known haploinsufficient genes (condition 6, not done). Combine this with the NA-observed fix rerun.
- **Run the ridge fit in the pipeline.** It's implemented on the local branch `jg/rr_ridge_fit` as the `ridge` feature-selection token. Locally it beat the current method with α, needed no α itself, and made the two DNM sets equal (see "Trying the alternatives" below). Next, run it on the cluster through `regional_rate_pipeline.py`.
- **Train regional effects on gnomAD rare variants.** The local version, trained on exome synonymous counts from the other chromosome half, scored best of all. The cleaner version trains on intronic rare variants from gnomAD genomes and needs the cluster.
- **Rebuild the background with gnomAD-observed sites included** (item 6). The approximate test moved the slope toward 1 but didn't lower the noise, so this is lower priority.
- **Recover the missing X/Y features**, or confirm they can't be recovered. Only this could make an X/Y correction worth more than about 1%.
- **Split the leftover noise into its two parts** (see "What the leftover noise is"). The site-level part belongs to the mutation model. The gene-level part could be reached by a nonlinear model trained on gnomAD rare variants, which removed 13% more locally. DNMs confirm most of the gene-level part is real mutation rate (see "What explains them besides size").
- **Report outliers on unflagged genes** (item 3 and the outlier figure).
- **Come back to the synonymous outliers.** Most of them have a real excess in their mutation rate, which DNMs confirm (see "What explains them besides size"). Open threads:
  - **Maternal hotspot regions:** the rr leaves about a third of their excess (8p: 1.31 to 1.20; DNMs 1.32). Try a Mb-scale hotspot feature, such as smoothed maternal C>G DNM density, or a regional term learned from gnomAD rare variants.
  - **Gene-level mutation rate:** DNMs confirm about two-thirds of the gene-to-gene spread the rr leaves. Consider whether a gene's own (shrunken) synonymous obs/exp should adjust its missense and LoF expected counts.
  - **The 11 "other" outliers** (HSP90AA1, EP300, MECP2, WFS1, MED16 and others): DNMs confirm the excess (1.73), but nothing measured explains it.
  - **Histone genes:** DNMs are 2.8 times expected. Check whether the extra variants and DNMs sit at positions where the paralogs differ, then decide whether to flag them.
  - **Repeat and segdup genes the QC flags miss** (FLG, MUC genes, olfactory receptors, PSG2, CT47B1, ACTG1, HSPA8): DNMs don't confirm their excess. Consider a flag.
  - **`low_exome_coverage` genes:** they all have too many variants, which suggests the coverage correction pulls expected counts too low. DNMs can't test this, because the DNM cohorts lose calls there too (DNM obs/exp 0.29).
  - **The outlier cutoff:** a size-aware rule, since a fixed |z| > 8 favours long genes.
- **To choose between the DNM sets:**
  - With a ridge fit, the two sets tie locally (2.475 vs 2.479), so the choice can rest on other grounds: more DNMs, chrX DNMs, and the newer annotation.
  - Run Ruchit's three diagnostics: downsample the new set to 17k, drop the 16 new features, and measure how much each feature carries over to gnomAD.
  - Look up the 282 dropped DNMs in the raw files, and try a union of the two sets if they were lost to re-annotation.
  - Compare the sets on unflagged-gene outliers and, after the full pipeline run, on LOEUF.
- **Mostly not worth chasing:** the remaining rise in LOEUF deciles 6–9. Binned by s_het, the top decile reaches only 1.04 instead of 1.11 (panel i), so most of the rise comes from how LOEUF sorts genes. Genes whose expected count is off in every variant class get pushed into the top LOEUF deciles. In unconstrained genes, synonymous and LoF obs/exp correlate (r = 0.22).

## Why more DNMs did worse

- Bonferroni is a fixed p-value cutoff, so 37% more DNMs let weaker features through. Of the 34 old features, 26 are kept and 8 dropped, and 16 new ones are added. Most of the new ones are GC, CpG-island and sperm-methylation features in non-CpG contexts.
- These features do predict DNMs: the held-out DNM gain rises from 86 to 132. But only about 20% of the change in rr shows up in gnomAD syn counts (slope 0.15–0.22, vs 0.67 for the old rr), and the loss is entirely in non-CpG contexts.
- At the stricter 1e-4 cutoff, which admits no new weak features, the new set performs the same as the old. With α fitted, the new set at Bonferroni comes close to the old set (a small gap remains). So the extra signal is real, but it is overstated at full strength.
- Why DNM-predictive features carry over poorly is untested. The candidates:
  - gnomAD expected counts already account for CpG methylation and coverage.
  - Polymorphism reflects long-term processes (GC-biased gene conversion, weak selection), not one generation of mutation.
  - How well exome DNM calling detects variants varies with GC differently from gnomAD's.
  - The background-selection effect in item 6.

## Recalibrating rr is normal

- **What the strength α is.** Raising rr to a power α multiplies every fitted feature effect by α. In statistics this is recalibration with a calibration slope, or shrinkage. It is not a fudge factor.
- **Why rr needs it.** There are two standard reasons a fitted model predicts differences that are too extreme, and rr has both:
  - **It's applied to data it wasn't trained on.** rr is fit on de novo mutations and applied to gnomAD polymorphism. A model moved to new data usually needs its calibration slope refit (Cox 1958; Steyerberg, *Clinical Prediction Models*).
  - **Its features are chosen by p-value, then fit without a penalty.** Effects of features picked that way are biased away from zero, so predictions need shrinking (Copas 1983; van Houwelingen & le Cessie 1990).
- **The rest of expected counts is already calibrated.** Mutation rates are estimated from genomes, then calibrated to exome synonymous sites by the plateau and coverage models. rr is the only part of the expected count used as estimated, without that calibration.
- **Not doing it before was probably an oversight.** Gnocchi multiplies expected counts by rr directly, with no calibration step, and the legacy `syn_adj_r` reused that code. Neither checked rr's calibration in coding sequence. Once checked, every version is too strong, legacy included.
- **The caveat.** A strength chosen on gnomAD synonymous counts will always look better on them. So choose it on held-out chromosomes, as Ruchit did, and judge it downstream (LOEUF after the full pipeline run).

## Trying the alternatives

![Out-of-fold noise and calibration slope for each alternative](followup/rr_alternatives.png)

The black row is the recommended ridge fit. In this figure it's fit on held-out chromosomes, so its numbers differ slightly from the all-chromosome version in the other figures.

How these were run:
- **Data:** Ruchit's training sites and evaluation export, run locally.
- **Held-out chromosomes:** every rr is fit on one chromosome half and scored on the other. That puts the numbers a bit above the all-chromosome ones elsewhere (his 17k rr: 2.60 here vs 2.56).

| Approach | DNMs | Noise as fitted | Calibration slope | Noise with best α | Syn outliers (no QC flag) |
|---|---|---|---|---|---|
| No correction | — | 2.908 | — | — | 270 (123) |
| Now: Bonferroni features, no penalty | 17k | 2.604 | 0.64 | 2.502 | 221 (71) |
| | Satterstrom | 2.759 | 0.55 | 2.563 | 232 (83) |
| Ridge on all features, penalty chosen on held-out DNMs | 17k | 2.475 | 0.80 | 2.484 | 218 (70) |
| | Satterstrom | 2.479 | 0.76 | 2.439 | 219 (70) |
| Sample split: pick features on half, fit on the other half | 17k | 2.577 | 1.06 | 2.593 | 247 (100) |
| | Satterstrom | 2.616 | 1.22 | 2.617 | 251 (102) |
| Background without gnomAD-observed sites* | 17k | 2.516 | 0.82 | 2.533 | 227 (76) |
| Background with all possible sites* | 17k | 2.534 | 0.93 | 2.596 | 233 (83) |
| No DNMs: fit on gnomAD synonymous counts | — | 2.365 | 0.86 | 2.377 | 216 (70) |

\*These two backgrounds are built from the exported MANE Select cells, so compare them only with each other. The Satterstrom set shows the same pattern: slope 0.72 → 0.89.

- **A penalty does α's job, and does it better.** Ridge on all features, with the penalty chosen from held-out DNMs alone (no gnomAD data), beats the current method with α added (2.48 vs 2.50) and doesn't need α itself.
- **The penalty also removes the cost of more DNMs.** With ridge the two sets tie (2.475 vs 2.479). The Satterstrom set only did worse because the fit had no penalty.
- **Picking and fitting features on the same data is why rr is too strong.** Picking them on one half and fitting on the other gives a slope of about 1 with no α. It's noisier because each step sees only half the data, so ridge gets the same benefit more cheaply.
- **The background design changes strength, not noise.** Including gnomAD-observed sites moves the slope toward 1 (0.82 → 0.93) but doesn't lower the noise in this approximation.
- **Fitting regional effects on gnomAD synonymous counts does best.** Trained on one chromosome half and scored on the other, it reaches 2.37 with the fewest outliers, about 4% below the best DNM-based rr. Two caveats:
  - It learns all regional variation in synonymous counts, not just mutation rate: coverage, mapping and GC-biased gene conversion too. That's arguably what expected counts need, but it assumes the same regional effects apply to LoF and missense.
  - Synonymous counts are both the training and the test data here, on different chromosomes. So the downstream LOEUF check matters most. A cleaner test would train on intronic rare variants from gnomAD genomes.

Code: `followup/rr_shrinkage_experiments.py`. Results: `followup/rr_shrinkage_experiments.tsv` and `.pass2.tsv`.

## Ruchit's 10/27 follow-up, and how ridge compares

His `rr_followup_tests_report.html` tests why more DNMs did worse:
- **It isn't the new cohorts.** Random 17k subsets of the new DNMs vary as much as old vs new does.
- **It's the features that pass Bonferroni.** With the old feature set, the new DNMs beat the old ones (slope 0.68 → 0.80).
- **The added features carry over to gnomAD poorly** (0.48), while the ones they replaced carried over well (1.13).
- **His conclusion: Bonferroni selection on DNMs is a noisy filter.** He still recommends the new DNMs with a 1e-2 cutoff at α 0.625, though he is unsure about adding a strength adjustment.
- **He reached two conclusions we had too:**
  - his core set selects on the outcome, since it drops outlier-flagged genes
  - the training background excludes gnomAD-observed sites

Ridge answers his conclusion directly: it removes the selection step, so no cutoff and no α are needed. Head to head, on his data and scored the same way (my code reproduces his numbers for both of his options exactly):

| Option | Noise, his core set | Noise, core set not selected on outcome | Outliers | Slope |
|---|---|---|---|---|
| No correction | 2.908 | 3.504 | 270 | — |
| Old DNMs, Bonferroni, α 1 (his reference) | 2.557 / 2.606 | 2.908 / 2.978 | 218 / 221 | 0.68 |
| New DNMs, 1e-2, α 0.625 (his recommendation) | 2.496 / 2.504 | 2.859 / 2.868 | 222 / 222 | 0.81 |
| **New DNMs, ridge, no α** | **2.459 / 2.476** | **2.775 / 2.789** | **218 / 219** | **0.86** |
| Old DNMs, ridge, no α | 2.458 / 2.490 | 2.789 / 2.826 | 223 / 222 | 0.93 |

Values are all-chromosome / out-of-fold models. The core set not selected on the outcome uses expected ≥ 20 with QC-flagged genes removed (17,706 genes). His version of it has 16,947 genes, so absolute values differ, but his differences between options match mine.

- **Ridge is lower on noise than his recommendation in every column,** by 0.03–0.04 on his core set and about 0.08 on the core set not selected on the outcome. It also has about the same or fewer outliers.
- **No cutoff and no α means two fewer settings to choose and defend**, which addresses his hesitation about the strength adjustment.
- **With ridge, the old and new DNM sets perform the same,** which matches his finding that the DNMs themselves aren't the problem.

Code: `followup/compare_ridge_vs_cutoff_alpha.py`. Results: `followup/ridge_vs_cutoff_alpha.tsv`.

## What the leftover noise is

![Leftover noise by gene size, and what gene-level covariates remove](followup/rr_leftover_noise.png)

This uses the 17,459 genes with no release flag of any kind and at least 20 expected synonymous variants; the flags already handle the rest. On this set, noise is 2.57 with no correction and 2.15 with the ridge rr, against at most 0.70 from chance.

- **About half of the leftover excess doesn't grow with gene size** (left). Even small genes sit well above chance (1.1–1.9). That part comes from individual sites, for example the three possible variants at a position sharing its coverage and mappability. No regional correction can remove it; a better per-site model could.
- **The other half grows with gene size**, as a gene-wide rate error of about 5% would. That is the part a regional or gene-level correction can reach.
- **It's concentrated in a few hundred genes.** The noisiest 5% of genes carry half the leftover excess, and 61% of the excess comes from genes with more variants than expected. The noisiest unflagged genes include:
  - histone cluster genes (synonymous o/e ≈ 1.8, despite clean coverage and mapping)
  - mucins and other tandem-repeat genes
  - HLA genes
  - a few highly constrained, broadly expressed genes (HNRNPA2B1, SFPQ, BRD2, USP34)
- **Gene-level covariates help a little** (right). Each model is fit on gnomAD synonymous counts from the other chromosome half.
  - Gene quality metrics, gene size and GTEx expression barely help on their own: at most 1.5%.
  - A nonlinear model of the same 13 regional features removes another 7%, so the per-context logistic rr isn't using everything these features carry.
  - Everything together removes 13% (2.15 → 1.87), about a fifth of what's left above chance.
- **Most of the gene-level part is real.** Synonymous DNMs follow gnomAD synonymous obs/exp across genes, at about two-thirds of the slope expected if all the spread above chance were mutation rate (see "What explains them besides size"). Histone genes are included: their DNMs are 2.8 times expected.

What this suggests:
1. Treat what's left as two problems. The site-level part belongs to the mutation model (per-site coverage, wider sequence context). The gene-level part belongs to the regional correction.
2. For the gene-level part, the most promising route is a model trained on gnomAD rare variants that uses the regional features flexibly (nonlinearly). This matches the "fit on gnomAD synonymous" result earlier.
3. Check histone clusters for paralog artifacts. Consider flagging repeat and segdup genes, whose excess DNMs don't confirm.
4. Expression breadth (a stand-in for transcription-coupled repair) adds little on its own. Replication timing is untested; it needs an external track.

Code: `followup/residual_noise.py`. Results: `followup/residual_noise.*.tsv`.

## Check run today: what the syn outliers are

![Synonymous outliers by release gene flag, and histone gene z-scores](followup/rr_outliers.png)

- **Left (obs ≪ exp):** under every synonymous rr, all but 0–1 of these carry `low_exome_mapping_quality`. They are segdup and paralog genes, with median syn o/e about 0.2. No rr can fix them, and none should.
- **Right (obs ≫ exp):** rr changes these, but not all of them.
  - **Histone genes (lavender):** 23–24 excess outliers under every version. None has a QC gene flag, though the release already marks all 24 `outlier_syn`. No rr removes any of them.
    - There are 72 histone genes in all: 54 in the chr6p22 cluster, 5 in the chr1q21 cluster, and 13 elsewhere. Every outlier is in the chr6p22 cluster, which is 24 of its 54 genes.
    - The excess covers the whole cluster, not just the outliers: pooled syn o/e is 1.88 in the chr6p22 cluster, 1.32 in chr1q21 and 1.23 elsewhere. 45 of the 72 genes have o/e above 1.5.
    - **Panel b:** the whole chr6p22 histone cluster is shifted to the right of every other gene: median z 7.4, against 0.1 for all other genes without a QC flag. The histone genes elsewhere sit in between, at a median of 2.3. The 24 outliers are just the part of that shifted distribution that crosses z = 8. The 33 non-histone genes in the same 25–28.5 Mb stretch sit slightly *below* expected (median z −1.1), the opposite of their histone neighbours.
    - **It's specific to histone genes, not the region.** The 33 non-histone genes interleaved in the same chr6 stretch (25–28.5 Mb) have pooled o/e 0.92, and the neighbouring MHC has 0.91. So neither regional mutation rate nor MHC diversity explains it, which is why no rr helps. Synonymous DNMs in the cluster are 2.8 times expected (see "What explains them besides size"), so a real mutation-rate excess is likely. Paralog mis-mapping could still affect both gnomAD and DNM calls.
  - **`low_exome_coverage` (light gray, right):** 18–20 genes with too many variants. It's the only release gene flag besides `low_exome_mapping_quality`. These genes have more variants than the coverage-corrected expected count.
  - **Everything else without a flag:** these drop from 98 to 42 with legacy `syn_adj_r`, 48–52 with Ruchit's tables, and 46 with the ridge rr.
- **Gnocchi `adj_r`** creates about 236 deficit outliers among genes without a flag, by pushing ordinary genes below their expected counts.
- **Unflagged genes only:** synonymous outliers fall from 123 to 66 (legacy), 71–76 (Ruchit's tables) and 71 (ridge), about half. After rr, histone genes are about a third of the unflagged outliers left.

- **Panels c and d, gene size:** for the same relative excess, z grows with the square root of a gene's expected count. Panel c plots expected count, which is what z depends on. Panel d plots coding length; there are about 160 expected synonymous variants per kb, so the two track each other closely.
  - Large genes cross |z| = 8 with only 20–50% more variants than expected. CSMD1, for example, is at obs/exp 1.49 with about 1,900 expected.
  - Histone genes have only about 60–130 expected, so they need nearly double.
  - The non-histone outliers are mostly large genes: median 599 expected, against 212 for all genes.
  - So outlier counts favour large genes, and a fixed |z| cutoff is not a fixed effect size.

**The 47 non-histone outliers left with the ridge rr** (46 with too many variants, plus MUC4 with too few). These genes carry no QC gene flag, but 45 carry a release outlier flag (41 `outlier_syn`, 37 `outlier_mis`), and 37 are outliers under all four synonymous rr versions.

### What explains them besides size

Size decides which genes cross |z| > 8, not why they have extra variants: gnomAD obs/exp doesn't change with gene size (1.01 to 1.00 across size fifths). For most of the 46, the mutation rate is genuinely higher than the model predicts.

![De novo mutations, maternal hotspot regions and ancestry checks for the outliers](followup/rr_outlier_drivers.png)

- **a. De novo mutations (DNMs) show the same excess.** DNMs measure the mutation rate directly, in different people, so they don't share gnomAD's population history or most of its artifacts. Here they are the Satterstrom 2026 synonymous DNMs, with obs/exp scaled so genes without a flag are at 1.

  | Outliers with too many variants | Genes | gnomAD obs/exp | DNM obs/exp (95% CI) |
  |---|---|---|---|
  | in a maternal hotspot region | 23 | 1.36 | 1.47 (1.21–1.73) |
  | other | 11 | 1.31 | 1.73 (1.33–2.16) |
  | repeat, segdup or low coverage | 12 | 1.41 | 1.15 (0.83–1.51) |
  | histone genes, chr6p22 cluster | 54 | 1.83 | 2.84 (2.16–3.56) |
  | genes of the same size without a flag | 1,808 | 1.00 | 1.05 (1.02–1.08) |

  - The excess holds in non-CpG contexts alone (DNMs 1.51), so it isn't from methylation, which the DNM expectation leaves out.
  - **Repeat genes are the exception:** DNMs don't confirm their excess. These are the ones likely to be mapping artifacts in gnomAD: FLG, MUC5AC, MUC3A, MUC17, the olfactory receptors, PSG2, CT47B1, ACTG1, HSPA8 and MROH1.
  - **Histone genes:** DNMs are 2.8 times expected (63 against 22), more than gnomAD shows. So the histone excess is not only a gnomAD paralog artifact.
- **b. Across all genes, DNMs follow synonymous obs/exp** (slope 0.51 ± 0.07). If all the gene-to-gene spread above chance were mutation rate, the slope would be about 0.75. So about two-thirds of the spread the rr leaves behind is real mutation-rate variation.
- **c. Half the outliers sit in known maternal mutation hotspots.** 23 of the 46 are in the regions where [Jónsson et al. 2017](https://doi.org/10.1038/nature24018) found raised maternal C>G DNM rates: the first 20 Mb of 8p, 2p16–22, 7p, 9p, 16p and 16q. Chance would put about 4 there (p = 3 × 10⁻¹³).
  - This explains the chr16 and chr8p clusters, and also adds ZFP36L2 and MSH6 (2p), DNAH11 (7p), and KCNV2 and KIAA2026 (9p).
  - [Goldmann et al. 2018](https://pubmed.ncbi.nlm.nih.gov/29507425/) link these regions to double-strand breaks in aging oocytes.
- **d. The rr only partly corrects these regions,** even for genes that aren't outliers.
  - Across the six regions, the ridge rr moves obs/exp only from 1.09 to 1.06. DNMs agree with what's left: 1.07 (1.02–1.12).
  - On 8p it moves from 1.31 to 1.20, and DNMs are at 1.32.
  - The rr has a maternal clustered-DNM feature, but it is zero for most of these genes.
- **e. The excess is the same in every ancestry group, and at every sample size from 1,000 to 500,000** (about 1.4–1.6 everywhere).
  - That rules out ancestry-specific causes such as copy-number or diversity differences.
  - It also rules out per-sample sequencing errors, which would grow with sample size.
  - It sits inside the DNM range for the same genes (gray band).
- **Other checks:**
  - **It covers the whole gene.** A median of 99% of each gene's 1-kb windows are above expected, against 50% in genes of the same size. Only FLG, MUC5AC, MUC3A, MROH1, WFS1, ARHGAP8, CT47B1 and MECP2 have it concentrated in part of the gene.
  - **It shows in every context class.** CpG sites show less (A/T 1.36, other C 1.46, CpG 1.11) because they are close to saturated in gnomAD.
  - **Missense rises by the same factor**, as a mutation-rate cause would predict.
  - **Compared with genes of the same size,** the outliers are closer to telomeres (median 9 against 23 Mb), and the rr already raises them slightly (1.02 against 0.98). Coverage, mapping quality, segdup fraction and testis expression don't differ.

Code: `followup/outlier_drivers.py`. Results: `followup/outlier_drivers.*.tsv`.

<details>
<summary>Table: the 47 genes, sorted by z</summary>

| Gene | Chr | Group | Syn obs | Syn exp (ridge) | Obs/exp | z (ridge) | DNMs obs / exp | Missense also high | Release outlier flag | Notes |
|---|---|---|---|---|---|---|---|---|---|---|
| CSMD1 | 8 | hotspot 8p | 2,862 | 1,919 | 1.49 | 21.5 | 13 / 8.8 | yes | outlier_syn, outlier_mis | 8p23–p22 |
| MUC5AC | 11 | repeat / segdup | 1,975 | 1,264 | 1.56 | 20.0 | 7 / 10.6 | yes | outlier_syn, outlier_mis | tandem repeats; segdup 0.57; AN90 0.12 |
| FLG | 1 | repeat / segdup | 2,551 | 1,829 | 1.39 | 16.9 | 7 / 7.2 | yes | outlier_syn, outlier_mis | tandem repeats; segdup 0.87 |
| PIEZO1 | 16 | hotspot 16q | 2,427 | 1,730 | 1.40 | 16.8 | 14 / 8.0 | yes | outlier_syn, outlier_mis | 16q23–24 |
| SRRM2 | 16 | hotspot 16p | 2,025 | 1,454 | 1.39 | 15.0 | 9 / 6.5 | yes | outlier_syn, outlier_mis | 16p13.3 |
| MUC3A | 7 | repeat / segdup | 1,398 | 954 | 1.46 | 14.4 | 8 / 5.3 | yes | outlier_syn, outlier_mis | tandem repeats; AN90 0.42 |
| MROH1 | 8 | repeat / segdup | 921 | 599 | 1.54 | 13.2 | 4 / 4.3 | yes | outlier_syn, outlier_mis | AN90 0.27 |
| DNAH11 | 7 | hotspot 7p | 2,461 | 1,908 | 1.29 | 12.7 | 8 / 7.4 | yes | outlier_syn, outlier_mis |  |
| GSE1 | 16 | hotspot 16q | 1,115 | 776 | 1.44 | 12.2 | 12 / 3.7 | yes | outlier_syn, outlier_mis | 16q23–24 |
| MAZ | 16 | hotspot 16p | 510 | 307 | 1.66 | 11.5 | 4 / 2.0 | no | outlier_syn | AN90 0.66 |
| CACNA1H | 16 | hotspot 16p | 2,202 | 1,730 | 1.27 | 11.3 | 14 / 8.6 | yes | outlier_syn, outlier_mis | 16p13.3 |
| FDFT1 | 8 | hotspot 8p | 351 | 199 | 1.77 | 10.8 | 1 / 0.8 | yes | outlier_syn, outlier_mis | 8p23–p22 |
| ZFP36L2 | 2 | hotspot 2p16–22 | 600 | 388 | 1.55 | 10.7 | 2 / 2.7 | yes | outlier_syn, outlier_mis |  |
| PCM1 | 8 | hotspot 8p | 1,085 | 787 | 1.38 | 10.6 | 0 / 2.8 | yes | outlier_syn, outlier_mis | 8p23–p22 |
| DNAH17 | 17 | other | 3,031 | 2,502 | 1.21 | 10.6 | 19 / 11.8 | no | outlier_syn |  |
| FANCA | 16 | hotspot 16q | 1,014 | 730 | 1.39 | 10.5 | 3 / 3.2 | yes | outlier_syn, outlier_mis | 16q23–24 |
| CTSB | 8 | hotspot 8p | 320 | 179 | 1.78 | 10.5 | 1 / 0.8 | yes | outlier_syn, outlier_mis | 8p23–p22 |
| RP1L1 | 8 | hotspot 8p | 1,813 | 1,419 | 1.28 | 10.5 | 10 / 6.0 | yes | outlier_syn, outlier_mis | 8p23–p22 |
| IRF2BPL | 14 | other | 768 | 529 | 1.45 | 10.4 | 2 / 3.4 | no | outlier_syn |  |
| WFS1 | 4 | other | 860 | 608 | 1.41 | 10.2 | 5 / 3.0 | yes | outlier_syn, outlier_mis |  |
| ACTG1 | 17 | repeat / segdup | 422 | 259 | 1.63 | 10.2 | 2 / 1.3 | no | outlier_syn | segdup 1.00 |
| SLC38A8 | 16 | hotspot 16q | 447 | 278 | 1.61 | 10.2 | 0 / 1.3 | yes | outlier_syn, outlier_mis | 16q23–24 |
| HSP90AA1 | 14 | other | 477 | 301 | 1.58 | 10.1 | 5 / 1.2 | no | outlier_syn |  |
| EP300 | 22 | other | 1,340 | 1,027 | 1.31 | 9.8 | 7 / 3.9 | no | outlier_syn |  |
| HSPA8 | 11 | repeat / segdup | 436 | 275 | 1.59 | 9.7 | 1 / 1.0 | no | outlier_syn | segdup 1.00 |
| DLC1 | 8 | hotspot 8p | 986 | 727 | 1.36 | 9.6 | 5 / 3.1 | yes | outlier_syn, outlier_mis | 8p23–p22 |
| OR4A5 | 11 | repeat / segdup | 245 | 134 | 1.83 | 9.6 | 3 / 0.5 | yes | outlier_syn, outlier_mis | paralog family; segdup 1.00 |
| MED16 | 19 | other | 897 | 653 | 1.37 | 9.5 | 8 / 3.2 | yes | outlier_syn, outlier_mis |  |
| ATP2C2 | 16 | hotspot 16q | 804 | 576 | 1.40 | 9.5 | 9 / 2.8 | yes | outlier_syn, outlier_mis | 16q23–24 |
| CTU2 | 16 | hotspot 16q | 482 | 316 | 1.53 | 9.3 | 1 / 1.5 | yes | outlier_syn, outlier_mis | 16q23–24 |
| MECP2 | X | other | 310 | 185 | 1.67 | 9.2 | 4 / 1.3 | no | outlier_syn |  |
| OR4C16 | 11 | repeat / segdup | 233 | 129 | 1.81 | 9.2 | 1 / 0.4 | yes | outlier_syn, outlier_mis | paralog family |
| WWOX | 16 | hotspot 16q | 359 | 223 | 1.61 | 9.1 | 2 / 1.1 | yes | outlier_syn, outlier_mis | 16q23–24 |
| PSG2 | 19 | repeat / segdup | 261 | 151 | 1.73 | 9.0 | 1 / 0.6 | yes | outlier_syn, outlier_mis | paralog family; segdup 1.00 |
| MUC17 | 7 | repeat / segdup | 2,325 | 1,931 | 1.20 | 9.0 | 10 / 7.2 | yes | outlier_mis | tandem repeats |
| PKD1L3 | 16 | hotspot 16q | 996 | 752 | 1.32 | 8.9 | 2 / 2.8 | yes | outlier_mis |  |
| DMD | X | other | 1,223 | 955 | 1.28 | 8.7 | 7 / 5.4 | yes | outlier_syn, outlier_mis |  |
| KCNV2 | 9 | hotspot 9p | 518 | 356 | 1.46 | 8.6 | 3 / 1.8 | yes | outlier_syn, outlier_mis |  |
| CARHSP1 | 16 | hotspot 16p | 161 | 83 | 1.95 | 8.6 | 0 / 0.4 | yes | outlier_syn, outlier_mis | 16p13.3 |
| CT47B1 | X | repeat / segdup | 228 | 133 | 1.72 | 8.3 | 0 / 1.0 | yes | outlier_syn, outlier_mis | paralog family; segdup 1.00 |
| ARHGAP8 | 22 | other | 375 | 245 | 1.53 | 8.3 | 2 / 1.2 | yes | outlier_syn, outlier_mis |  |
| OR4C46 | 11 | repeat / segdup | 224 | 130 | 1.72 | 8.2 | 2 / 0.4 | yes | outlier_mis | paralog family; segdup 1.00 |
| MSH6 | 2 | hotspot 2p16–22 | 830 | 625 | 1.33 | 8.2 | 2 / 2.4 | yes | outlier_syn, outlier_mis |  |
| KIAA2026 | 9 | hotspot 9p | 1,222 | 967 | 1.26 | 8.2 | 5 / 3.8 | yes | outlier_syn, outlier_mis |  |
| LRIG1 | 3 | other | 772 | 577 | 1.34 | 8.1 | 5 / 2.7 | no | — |  |
| ID1 | 20 | other | 176 | 97 | 1.82 | 8.1 | 1 / 0.6 | yes | outlier_mis |  |
| MUC4 | 3 | repeat / segdup | 2,351 | 2,874 | 0.82 | -9.8 | 23 / 13.6 | no | — | tandem repeats; segdup 0.63; AN90 0.69 |

"Missense also high" means release missense z < −8, which is the release's sign for too many variants. "DNMs obs / exp" counts Satterstrom 2026 synonymous DNMs in the gene against the number expected from context rates and the ridge rr; per-gene counts are small, so read them in groups. Group: in a maternal hotspot region; repeat or segdup (segdup fraction ≥ 0.5, AN90 < 0.5, or a MUC, OR or PSG gene); or other. Data: `followup/ridge_unflagged_outliers.tsv`, `followup/outlier_drivers.genes.tsv`.

</details>

Histone genes are assigned by gene symbol (H1-, H2AC, H2BC, H3C, H4C, H2AW, H2BW, H3- and H4-). A histone transcript with a release flag is counted under its flag.

The counts are recomputed from Ruchit's evaluation export: all 18,406 Ensembl MANE Select transcripts including X/Y, raw |z| > 8, and his calibrated tables. Flags come from the v4.1.2 release. The totals reproduce his reported counts (270, 447, 217, 219 and 224). Code: `followup/remap_rr_to_shet.py --outliers` (histone pattern in `GENE_FAMILIES`).
