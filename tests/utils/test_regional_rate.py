"""Tests for the regional mutation rate (rr) functions in ``gnomad_constraint.utils.regional_rate``."""

from typing import Optional

import hail as hl
import numpy as np
import pytest
from gnomad.utils.constraint import MODELS_HT_SCHEMA

from gnomad_constraint.resources.constants import (
    GENOMIC_FEATURE_COLUMNS,
    RR_CONTEXT_FEATURES,
    RR_CONTEXTS,
)
from gnomad_constraint.utils.regional_rate import (
    compute_rr_ht,
    create_rr_site_features_ht,
    get_rr_fold_expr,
    get_rr_normalization_sites_ht,
    remove_qc_flagged_genes,
    select_rr_features,
)

INTERVAL_TYPE = hl.tinterval(hl.tlocus("GRCh38"))


def _locus(pos: int, contig: str = "chr1") -> hl.Locus:
    """Return a GRCh38 locus."""
    return hl.Locus(contig, pos, "GRCh38")


class TestRrContextFeatures:
    """Test the RR_CONTEXT_FEATURES constant."""

    def test_cpg_contexts_leave_out_methylation_features(self) -> None:
        """
        Test that only CpG contexts leave out the methylation-correlated features.

        Every context has all 13 features at 4 scales (52), except the 4 CpG
        contexts, which leave out 5 features at every scale (32).
        """
        assert len(RR_CONTEXT_FEATURES) == 32
        counts = {c: len(f) for c, f in RR_CONTEXT_FEATURES.items()}
        cpg = {"ACG", "CCG", "GCG", "TCG"}
        assert {c for c, n in counts.items() if n == 32} == cpg
        assert all(n == 52 for c, n in counts.items() if c not in cpg)
        assert RR_CONTEXT_FEATURES["AAA"] == GENOMIC_FEATURE_COLUMNS
        assert "GC_content_1k" not in RR_CONTEXT_FEATURES["TCG"]
        assert "GC_content_1k" in RR_CONTEXT_FEATURES["CCA"]


def _metrics_ht(flagged: bool = True, span: tuple = (150, 250)) -> hl.Table:
    """
    Return constraint metrics rows with gene flags, all on chr1.

    ENST1 (positions ``span``) is flagged only if ``flagged``. ENST2 isn't flagged,
    and the flagged NM_1, non-MANE ENST3, and span-less ENST4 don't count.
    """
    flags = {"low_exome_coverage"} if flagged else set()
    metrics = [
        ("ENST1", True, flags, *span),
        ("ENST2", True, set(), 50, 120),
        ("NM_1", True, {"low_exome_coverage"}, 280, 320),
        ("ENST3", False, {"low_exome_coverage"}, 90, 110),
        ("ENST4", True, {"low_exome_coverage"}, None, None),
    ]
    return hl.Table.parallelize(
        [
            {
                "transcript": transcript,
                "mane_select": mane_select,
                "gene_flags": flags,
                "chromosome": "chr1",
                "start_position": start,
                "end_position": end,
            }
            for transcript, mane_select, flags, start, end in metrics
        ],
        hl.tstruct(
            transcript=hl.tstr,
            mane_select=hl.tbool,
            gene_flags=hl.tset(hl.tstr),
            chromosome=hl.tstr,
            start_position=hl.tint32,
            end_position=hl.tint32,
        ),
    )


class TestCreateRrSiteFeaturesHt:
    """Test the create_rr_site_features_ht function."""

    PREPROCESSED_TYPE = hl.tstruct(
        locus=hl.tlocus("GRCh38"),
        alleles=hl.tarray(hl.tstr),
        context=hl.tstr,
        vep=hl.tstruct(most_severe_consequence=hl.tstr),
        filters=hl.tstruct(exomes=hl.tset(hl.tstr), genomes=hl.tset(hl.tstr)),
        exomes_coverage=hl.tint32,
    )
    """Fields of the preprocessed Table the function reads."""

    @staticmethod
    def _globals(
        upper_cov_cutoff: Optional[int] = None,
    ) -> hl.expr.StructExpression:
        """Return preprocessed globals with the plateau models' coverage range."""
        return hl.struct(
            version="test",
            build_models_globals=hl.struct(
                high_cov_cutoff=90,
                upper_cov_cutoff=hl.or_missing(
                    upper_cov_cutoff is not None, hl.int32(upper_cov_cutoff or 0)
                ),
            ),
            apply_models_globals=hl.struct(exome_coverage_metric="AN_percent"),
        )

    def _preprocessed_ht(self, upper_cov_cutoff: Optional[int] = None) -> hl.Table:
        """
        Return preprocessed Table rows at chr1 positions 10-110.

        10 and 20 are synonymous DNMs, 20 also in exomes; 10 has a second
        synonymous allele absent from exomes. 30 has two synonymous alleles absent
        from exomes, 40 is in exomes and passes, 50 is in exomes as AC0, 60 is
        missense, and 70 is absent from exomes. The sites have different contexts,
        so a context mix-up shows. 80 (a DNM) and 90 are below the high coverage
        cutoff of 90, 100 has no coverage, and 110 has coverage 100.
        """
        rows = [
            (10, ["A", "G"], "AAA", "synonymous_variant", None, 95),
            (10, ["A", "C"], "AAA", "synonymous_variant", None, 95),
            (20, ["C", "T"], "ACG", "synonymous_variant", set(), 90),
            (30, ["A", "G"], "CAT", "synonymous_variant", None, 95),
            (30, ["A", "C"], "CAT", "synonymous_variant", None, 95),
            (40, ["A", "G"], "AAA", "synonymous_variant", set(), 95),
            (50, ["A", "G"], "AAA", "synonymous_variant", {"AC0"}, 95),
            (60, ["C", "T"], "ACG", "missense_variant", None, 95),
            (70, ["A", "G"], "TAC", "synonymous_variant", None, 95),
            (80, ["A", "G"], "AAA", "synonymous_variant", None, 50),
            (90, ["A", "G"], "AAA", "synonymous_variant", None, 89),
            (100, ["A", "G"], "AAA", "synonymous_variant", None, None),
            (110, ["A", "G"], "AAA", "synonymous_variant", None, 100),
        ]
        return hl.Table.parallelize(
            [
                {
                    "locus": _locus(pos),
                    "alleles": alleles,
                    "context": context,
                    "vep": {"most_severe_consequence": csq},
                    "filters": {"exomes": filters, "genomes": None},
                    "exomes_coverage": coverage,
                }
                for pos, alleles, context, csq, filters, coverage in rows
            ],
            self.PREPROCESSED_TYPE,
            key=["locus", "alleles"],
            globals=self._globals(upper_cov_cutoff),
        )

    @staticmethod
    def _dnm_ht() -> hl.Table:
        """Return DNMs at positions 10, 20, 60, and 80, and an indel at 70."""
        return hl.Table.parallelize(
            [
                {"locus": _locus(pos), "alleles": alleles}
                for pos, alleles in [
                    (10, ["A", "G"]),
                    (20, ["C", "T"]),
                    (60, ["C", "T"]),
                    (70, ["A", "AT"]),
                    (80, ["A", "G"]),
                ]
            ],
            hl.tstruct(locus=hl.tlocus("GRCh38"), alleles=hl.tarray(hl.tstr)),
            key=["locus", "alleles"],
        )

    @staticmethod
    def _features_ht() -> hl.Table:
        """Return two windows, positions 1-25 (features 1.0) and 26-1000 (2.0)."""
        return hl.Table.parallelize(
            [
                {
                    "interval": hl.Interval(
                        _locus(start), _locus(end), includes_end=True
                    ),
                    **dict.fromkeys(GENOMIC_FEATURE_COLUMNS, value),
                }
                for start, end, value in [(1, 25, 1.0), (26, 1000, 2.0)]
            ],
            hl.tstruct(
                interval=INTERVAL_TYPE,
                **dict.fromkeys(GENOMIC_FEATURE_COLUMNS, hl.tfloat64),
            ),
            key="interval",
        )

    def _site_features_ht(
        self,
        background_fraction: float,
        seed: int = 0,
        gene_flags_ht: Optional[hl.Table] = None,
        upper_cov_cutoff: Optional[int] = None,
    ) -> hl.Table:
        """Return the site features Table of the test rows."""
        return create_rr_site_features_ht(
            self._preprocessed_ht(upper_cov_cutoff),
            self._dnm_ht(),
            self._features_ht(),
            background_fraction,
            seed,
            gene_flags_ht=gene_flags_ht,
            gene_flags_version="4.1.2",
        )

    def test_sites(self) -> None:
        """
        Test the DNM and background sites.

        DNMs are the synonymous context sites in the DNM Table, even if in exomes
        (20); the missense DNM (60) and the indel (70) aren't. The background (all of
        it, with a fraction of 1) is the synonymous sites with no exomes filters at
        the other loci, one per locus (30, 70): sites in exomes, passing (40) or AC0
        (50), aren't in it, and neither is the second allele at the DNM locus 10.
        Sites below the plateau models' high coverage cutoff (the DNM at 80 and 90)
        or with no coverage (100) aren't in either set. Every site gets its window's
        features, and the preprocessed Table's globals are replaced by the record of
        how the sites were chosen.
        """
        ht = self._site_features_ht(1.0)
        assert hl.eval(ht.globals) == hl.Struct(
            rr_training_sites=hl.Struct(
                background_fraction=1.0,
                background_seed=0,
                qc_flagged_gene_filter=False,
                gene_flags_version=None,
                exome_coverage_metric="AN_percent",
                high_cov_cutoff=90,
                upper_cov_cutoff=None,
            )
        )
        assert list(ht.key) == ["locus"]
        rows = ht.collect()
        assert sorted(
            (r.locus.position, r.alleles, r.context, r.y, r.dist2telo_1k) for r in rows
        ) == [
            (10, ["A", "G"], "AAA", 1, 1.0),
            (20, ["C", "T"], "ACG", 1, 1.0),
            (30, ["A", "C"], "CAT", 0, 2.0),
            (70, ["A", "G"], "TAC", 0, 2.0),
            (110, ["A", "G"], "AAA", 0, 2.0),
        ]

    def test_upper_coverage_cutoff(self) -> None:
        """
        Test that sites above the plateau models' upper coverage cutoff are left out.

        With an upper cutoff of 99, the background site at 110 (coverage 100) goes,
        and the cutoff is recorded.
        """
        ht = self._site_features_ht(1.0, upper_cov_cutoff=99)
        assert hl.eval(ht.rr_training_sites.upper_cov_cutoff) == 99
        assert sorted(r.locus.position for r in ht.collect()) == [10, 20, 30, 70]

    def test_no_coverage_metric_raises(self) -> None:
        """
        Test that a preprocessed Table without its exome coverage metric raises.

        Guards against a preprocessed Table written before ``apply_models_globals``
        had ``exome_coverage_metric`` failing with an AttributeError only after the
        scan, as it did on a ``--test`` run.
        """
        preprocessed_ht = self._preprocessed_ht().annotate_globals(
            apply_models_globals=hl.struct(skip_coverage_model=False)
        )
        with pytest.raises(ValueError, match="Rerun prepare_data.py --preprocess"):
            create_rr_site_features_ht(
                preprocessed_ht, self._dnm_ht(), self._features_ht(), 1.0, 0
            )

    def test_qc_flagged_genes_left_out(self) -> None:
        """
        Test that with gene flags, the sites in a flagged gene's span are left out.

        The flagged gene spans 15-35, so the DNM at 20 and the background site at 30
        both go.
        """
        ht = self._site_features_ht(1.0, gene_flags_ht=_metrics_ht(span=(15, 35)))
        rr_training_sites = hl.eval(ht.rr_training_sites)
        assert rr_training_sites.qc_flagged_gene_filter
        assert rr_training_sites.gene_flags_version == "4.1.2"
        rows = ht.collect()
        assert sorted((r.locus.position, r.y) for r in rows) == [
            (10, 1),
            (70, 0),
            (110, 0),
        ]

    def test_no_background(self) -> None:
        """Test that a background fraction of 0 keeps only the DNM sites."""
        rows = self._site_features_ht(0.0).collect()
        assert sorted((r.locus.position, r.y) for r in rows) == [(10, 1), (20, 1)]

    def test_background_sample(self) -> None:
        """
        Test that the background is sampled per locus, at the fraction, by the seed.

        Each of 60 loci has 3 synonymous alleles absent from exomes. Sampling half
        the loci keeps about 30; sampling alleles before keeping one per locus would
        keep about 52. Different seeds give different samples, and a seed repeats.
        """
        n = 60
        context_ht = hl.Table.parallelize(
            [
                {
                    "locus": _locus(10 * (i + 1)),
                    "alleles": ["A", alt],
                    "context": "AAA",
                    "vep": {"most_severe_consequence": "synonymous_variant"},
                    "filters": {"exomes": None, "genomes": None},
                    "exomes_coverage": 95,
                }
                for i in range(n)
                for alt in ["C", "G", "T"]
            ],
            self.PREPROCESSED_TYPE,
            key=["locus", "alleles"],
            globals=self._globals(),
        )
        dnm_ht = self._dnm_ht().filter(False)

        def _background(fraction: float, seed: int) -> list:
            ht = create_rr_site_features_ht(
                context_ht, dnm_ht, self._features_ht(), fraction, seed
            )
            return sorted(r.locus.position for r in ht.collect())

        sample = _background(0.5, 0)
        assert len(set(sample)) == len(sample)
        assert 18 <= len(sample) <= 42
        assert sample == _background(0.5, 0)
        assert sample != _background(0.5, 1)
        assert len(_background(1.0, 0)) == n


class TestRemoveQcFlaggedGenes:
    """Test the remove_qc_flagged_genes function."""

    @pytest.mark.parametrize("flagged", [True, False])
    def test_flagged_ensembl_mane_select_spans(self, flagged: bool) -> None:
        """
        Test that only the spans of flagged Ensembl MANE Select transcripts are removed.

        A span includes its start and end. A transcript with no span is skipped. With
        no flagged gene, every site is kept.
        """
        positions = [100, 149, 150, 200, 250, 251, 300]
        ht = hl.Table.parallelize(
            [{"locus": _locus(p)} for p in positions],
            hl.tstruct(locus=hl.tlocus("GRCh38")),
            key="locus",
        )
        metrics_ht = _metrics_ht(flagged)
        out = remove_qc_flagged_genes(ht, metrics_ht).collect()
        expected = [100, 149, 251, 300] if flagged else positions
        assert [r.locus.position for r in out] == expected


class TestGetRrFoldExpr:
    """Test the get_rr_fold_expr function."""

    def test_folds(self) -> None:
        """
        Test that odd-numbered chromosomes are in the odd fold and even ones even.

        chrX gets a fold too, but has no features, so it never reaches a fit.
        """
        contigs = ["chr1", "chr2", "chr3", "chr20", "chrX"]
        ht = hl.Table.parallelize(
            [{"locus": _locus(10, c)} for c in contigs],
            hl.tstruct(locus=hl.tlocus("GRCh38")),
        )
        ht = ht.annotate(fold=get_rr_fold_expr(ht.locus))
        assert {r.locus.contig: r.fold for r in ht.collect()} == {
            "chr1": "odd",
            "chr2": "even",
            "chr3": "odd",
            "chr20": "even",
            "chrX": "odd",
        }


def _per_variant_ht(rows: list) -> hl.Table:
    """
    Return a per-variant expected Table, one row per site and transcript.

    Each row is (position, transcript, mane_select, annotation, context, adj
    expected count); the adj counts are at index 1 of ``expected_variants``.
    """
    ht = hl.Table.parallelize(
        [
            {
                "locus": _locus(pos),
                "alleles": ["A", "G"],
                "transcript": transcript,
                "mane_select": mane_select,
                "annotation": annotation,
                "context": context,
                "expected_variants": [100.0, exp],
            }
            for pos, transcript, mane_select, annotation, context, exp in rows
        ],
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            transcript=hl.tstr,
            mane_select=hl.tbool,
            annotation=hl.tstr,
            context=hl.tstr,
            expected_variants=hl.tarray(hl.tfloat64),
        ),
        key=["locus", "alleles"],
    )
    return ht.annotate_globals(
        exomes_freq_meta=[{"group": "raw"}, {"group": "adj"}],
        apply_models_globals=hl.struct(version="test"),
    )


class TestSelectRrFeatures:
    """Test the select_rr_features function."""

    @staticmethod
    def _screen_ht() -> hl.Table:
        """
        Return a feature screen of AAA, ACG (a CpG context), and CCC.

        A missing p-value is a screen fit that didn't converge.
        """
        rows = [
            ("AAA", "GC_content_1k", 1e-3),
            ("AAA", "GC_content_10k", 5e-4),
            ("AAA", "GC_content_100k", 1e-6),
            ("AAA", "LCR_1M", 9e-4),
            ("AAA", "SINE_1k", None),
            ("AAA", "SINE_10k", 2e-3),
            ("ACG", "LCR_1k", 1.2e-3),
            ("ACG", "GC_content_1k", 1e-9),
            ("CCC", "LINE_1k", 0.5),
        ]
        return hl.Table.parallelize(
            [{"context": c, "feature": f, "p_value": p} for c, f, p in rows],
            hl.tstruct(context=hl.tstr, feature=hl.tstr, p_value=hl.tfloat64),
            key=["context", "feature"],
        )

    def test_bonferroni(self) -> None:
        """
        Test Gnocchi's Bonferroni selection: p <= 0.05 / 4 / 13, or / 8 in CpG contexts.

        Only the smallest selected scale of a feature is kept (GC_content at 10k,
        not 100k; 1k is just above the cutoff). ACG's LCR passes its CpG cutoff,
        1.56e-3, which the non-CpG cutoff, 9.6e-4, would not. GC_content isn't one of
        a CpG context's features, so ACG doesn't get it whatever its p-value, a fit
        with no p-value is never selected, and CCC, with nothing selected, is left
        out.
        """
        assert select_rr_features(self._screen_ht(), "bonferroni") == {
            "AAA": ["GC_content_10k", "LCR_1M"],
            "ACG": ["LCR_1k"],
        }

    def test_p_value_cutoff(self) -> None:
        """Test that any other value is the p-value cutoff itself, for every context."""
        assert select_rr_features(self._screen_ht(), "1e-2") == {
            "AAA": ["GC_content_1k", "SINE_10k", "LCR_1M"],
            "ACG": ["LCR_1k"],
        }


class TestGetRrNormalizationSitesHt:
    """Test the get_rr_normalization_sites_ht function."""

    def test_synonymous_ensembl_mane_select_sites(self) -> None:
        """
        Test that only synonymous sites in Ensembl MANE Select transcripts are kept.

        The RefSeq MANE Select row of a site would count it twice. Each site keeps
        its adj expected count, found through ``exomes_freq_meta``.
        """
        ht = get_rr_normalization_sites_ht(
            _per_variant_ht(
                [
                    (10, "ENST1", True, "synonymous_variant", "AAA", 2.0),
                    (10, "NM_1", True, "synonymous_variant", "AAA", 2.0),
                    (20, "ENST2", False, "synonymous_variant", "ACG", 1.0),
                    (30, "ENST3", True, "missense_variant", "AAA", 5.0),
                    (40, "ENST4", True, "synonymous_variant", "ACG", 3.0),
                ]
            )
        )
        assert list(ht.row) == ["locus", "alleles", "context", "exp"]
        assert [(r.locus.position, r.context, r.exp) for r in ht.collect()] == [
            (10, "AAA", 2.0),
            (40, "ACG", 3.0),
        ]


class TestComputeRrHt:
    """Test the compute_rr_ht function."""

    @staticmethod
    def _features_ht() -> hl.Table:
        """Return three windows; the third has no feature value."""
        return hl.Table.parallelize(
            [
                {"interval": hl.Interval(_locus(s), _locus(e)), "a": a}
                for s, e, a in [(1, 1001, 0.5), (1001, 2001, 1.5), (2001, 3001, None)]
            ],
            hl.tstruct(interval=INTERVAL_TYPE, a=hl.tfloat64),
            key="interval",
        )

    @staticmethod
    def _models_ht(
        context_feature: Optional[dict] = None, converged: bool = True
    ) -> hl.Table:
        """
        Return models with rate ``expit(x - 1) / 0.5`` of one feature, penalty 0.3.

        :param context_feature: Each model's context and feature. Default is one AAA
            model of ``a``.
        :param converged: Whether the models converged. Default is True.
        """
        ht = hl.Table.parallelize(
            [
                {
                    "context": context,
                    "features": [feature],
                    "mean": [1.0],
                    "sd": [1.0],
                    "coefs": [1.0],
                    "intercept": 0.0,
                    "ave": 0.5,
                    "n_mutated": 10,
                    "n_background": 10,
                    "converged": converged,
                }
                for context, feature in (context_feature or {"AAA": "a"}).items()
            ],
            MODELS_HT_SCHEMA,
            key="context",
        )
        return ht.annotate_globals(
            ridge_penalty=0.3,
            rr_training_sites=hl.struct(background_seed=0),
            rr_feature_selection="all",
        )

    def test_rr(self) -> None:
        """
        Test rr per window and context, its normalization, and its globals.

        The window without a feature value has no rr and is dropped, so the site in
        it is left out of the normalization. AAA's rr, weighted by the sites' adj
        expected counts, averages 1 after normalization, and the unmodeled contexts
        stay 1.0.
        """
        ht = compute_rr_ht(
            self._features_ht(),
            self._models_ht(),
            _per_variant_ht(
                [
                    (10, "ENST1", True, "synonymous_variant", "AAA", 1.0),
                    (20, "ENST2", True, "synonymous_variant", "ACG", 2.0),
                    (1500, "ENST3", True, "synonymous_variant", "AAA", 3.0),
                    (2500, "ENST4", True, "synonymous_variant", "AAA", 5.0),
                ]
            ),
        )
        rr = {r.interval.start.position: r.adj_r for r in ht.collect()}
        assert sorted(rr) == [1, 1001]
        assert set(rr[1]) == set(RR_CONTEXTS)
        raw = {1: 2 / (1 + np.exp(0.5)), 1001: 2 / (1 + np.exp(-0.5))}
        factor = (1 * raw[1] + 3 * raw[1001]) / 4
        assert rr[1]["AAA"] == pytest.approx(raw[1] / factor)
        assert rr[1001]["AAA"] == pytest.approx(raw[1001] / factor)
        assert (1 * rr[1]["AAA"] + 3 * rr[1001]["AAA"]) / 4 == pytest.approx(1.0)
        assert rr[1]["ACG"] == 1.0
        rr_globals = hl.eval(ht.rr_globals)
        assert rr_globals.ridge_penalty == 0.3
        assert rr_globals.feature_selection == "all"
        assert rr_globals.normalized
        assert rr_globals.normalization_factors == pytest.approx(
            {"AAA": factor, "ACG": 1.0}
        )
        assert rr_globals.per_variant_apply_models_globals == hl.Struct(version="test")
        assert rr_globals.per_variant_row_count == 4
        assert rr_globals.training_sites == hl.Struct(background_seed=0)
        assert rr_globals.non_converged_contexts == []

    def test_without_normalization(self) -> None:
        """
        Test that without ``normalize``, rr is the models' rate and the globals say so.

        This is how Gnocchi's rr was used. No normalization factors are computed.
        """
        ht = compute_rr_ht(
            self._features_ht(),
            self._models_ht(),
            _per_variant_ht([(10, "ENST1", True, "synonymous_variant", "AAA", 1.0)]),
            normalize=False,
        )
        rr = {r.interval.start.position: r.adj_r for r in ht.collect()}
        assert rr[1]["AAA"] == pytest.approx(2 / (1 + np.exp(0.5)))
        assert rr[1001]["AAA"] == pytest.approx(2 / (1 + np.exp(-0.5)))
        rr_globals = hl.eval(ht.rr_globals)
        assert not rr_globals.normalized
        assert rr_globals.normalization_factors == {}
        assert rr_globals.per_variant_row_count == 1

    def test_non_converged_model_left_out(self) -> None:
        """
        Test that a context whose model didn't converge gets rr 1.0, and is recorded.

        Regression test: a non-converged fit was applied as is, and nothing past the
        models Table recorded which contexts it was.
        """
        ht = compute_rr_ht(
            self._features_ht(),
            self._models_ht(converged=False),
            _per_variant_ht(
                [(10, "ENST1", True, "synonymous_variant", "AAA", 1.0)],
            ),
        )
        assert all(r.adj_r["AAA"] == 1.0 for r in ht.collect())
        assert hl.eval(ht.rr_globals.non_converged_contexts) == ["AAA"]

    def test_window_missing_a_models_feature(self) -> None:
        """
        Test that a window missing one model's feature has no rr in that context.

        The second window has ``a`` but not ``b``, so it has an rr in AAA but not in
        ACA, whose model uses ``b``. Its ACA rr stays missing after normalization, so
        its expected counts in ACA aren't adjusted, and its ACA site is left out of
        ACA's factor. It used to get the rate of the average window, 1.0, which
        entered the factor and was then divided by it.
        """
        features_ht = hl.Table.parallelize(
            [
                {"interval": hl.Interval(_locus(s), _locus(e)), "a": a, "b": b}
                for s, e, a, b in [(1, 1001, 0.5, 2.0), (1001, 2001, 1.5, None)]
            ],
            hl.tstruct(interval=INTERVAL_TYPE, a=hl.tfloat64, b=hl.tfloat64),
            key="interval",
        )
        ht = compute_rr_ht(
            features_ht,
            self._models_ht({"AAA": "a", "ACA": "b"}),
            _per_variant_ht(
                [
                    (10, "ENST1", True, "synonymous_variant", "ACA", 1.0),
                    (1500, "ENST2", True, "synonymous_variant", "ACA", 3.0),
                ]
            ),
        )
        rr = {r.interval.start.position: r.adj_r for r in ht.collect()}
        assert sorted(rr) == [1, 1001]
        assert rr[1001]["AAA"] == pytest.approx(2 / (1 + np.exp(-0.5)))
        assert rr[1001]["ACA"] is None
        assert rr[1]["ACA"] == pytest.approx(1.0)
        assert hl.eval(ht.rr_globals.normalization_factors)["ACA"] == pytest.approx(
            2 / (1 + np.exp(-1.0))
        )
