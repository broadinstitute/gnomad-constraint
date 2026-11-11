"""Tests for the regional mutation rate (rr) functions in ``gnomad_constraint.utils.regional_rate``."""

import hail as hl
import pytest

from gnomad_constraint.resources.constants import (
    GENOMIC_FEATURE_COLUMNS,
    RR_CONTEXT_FEATURES,
)
from gnomad_constraint.utils.regional_rate import (
    create_rr_training_ht,
    get_rr_fold_expr,
    remove_qc_flagged_genes,
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


class TestCreateRrTrainingHt:
    """Test the create_rr_training_ht function."""

    @staticmethod
    def _context_ht() -> hl.Table:
        """
        Return context Table rows at chr1 positions 10-70.

        10 and 20 are synonymous DNMs, 20 also in exomes; 10 has a second
        synonymous allele absent from exomes. 30 has two synonymous alleles absent
        from exomes, 40 is in exomes and passes, 50 is in exomes as AC0, 60 is
        missense, and 70 is absent from exomes.
        """
        rows = [
            (10, ["A", "G"], "synonymous_variant", None),
            (10, ["A", "C"], "synonymous_variant", None),
            (20, ["C", "T"], "synonymous_variant", set()),
            (30, ["A", "G"], "synonymous_variant", None),
            (30, ["A", "C"], "synonymous_variant", None),
            (40, ["A", "G"], "synonymous_variant", set()),
            (50, ["A", "G"], "synonymous_variant", {"AC0"}),
            (60, ["C", "T"], "missense_variant", None),
            (70, ["A", "G"], "synonymous_variant", None),
        ]
        return hl.Table.parallelize(
            [
                {
                    "locus": _locus(pos),
                    "alleles": alleles,
                    "context": "AAA",
                    "vep": {"most_severe_consequence": csq},
                    "filters": {"exomes": filters, "genomes": None},
                }
                for pos, alleles, csq, filters in rows
            ],
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                context=hl.tstr,
                vep=hl.tstruct(most_severe_consequence=hl.tstr),
                filters=hl.tstruct(exomes=hl.tset(hl.tstr), genomes=hl.tset(hl.tstr)),
            ),
            key=["locus", "alleles"],
            globals=hl.Struct(version="test"),
        )

    @staticmethod
    def _dnm_ht() -> hl.Table:
        """Return DNMs at positions 10, 20, and 60, and an indel at 70."""
        return hl.Table.parallelize(
            [
                {"locus": _locus(pos), "alleles": alleles}
                for pos, alleles in [
                    (10, ["A", "G"]),
                    (20, ["C", "T"]),
                    (60, ["C", "T"]),
                    (70, ["A", "AT"]),
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

    def _training_ht(self, background_fraction: float, seed: int = 0) -> hl.Table:
        """Return the training Table of the test rows."""
        return create_rr_training_ht(
            self._context_ht(),
            self._dnm_ht(),
            self._features_ht(),
            background_fraction,
            seed,
        )

    def test_training_sites(self) -> None:
        """
        Test the DNM and background sites of the training set.

        DNMs are the synonymous context sites in the DNM Table, even if in exomes
        (20); the missense DNM (60) and the indel (70) aren't. The background (all of
        it, with a fraction of 1) is the synonymous sites with no exomes filters at
        the other loci, one per locus (30, 70): sites in exomes, passing (40) or AC0
        (50), aren't in it, and neither is the second allele at the DNM locus 10.
        Every site gets its window's features, and the context Table's globals are
        dropped.
        """
        ht = self._training_ht(1.0)
        assert len(ht.globals.dtype) == 0
        rows = ht.collect()
        assert sorted((r.locus.position, r.y, r.dist2telo_1k) for r in rows) == [
            (10, 1, 1.0),
            (20, 1, 1.0),
            (30, 0, 2.0),
            (70, 0, 2.0),
        ]

    def test_no_background(self) -> None:
        """Test that a background fraction of 0 keeps only the DNM sites."""
        rows = self._training_ht(0.0).collect()
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
                }
                for i in range(n)
                for alt in ["C", "G", "T"]
            ],
            self._context_ht().row.dtype,
            key=["locus", "alleles"],
        )
        dnm_ht = self._dnm_ht().filter(False)

        def _background(fraction: float, seed: int) -> list:
            ht = create_rr_training_ht(
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
        flags = {"low_exome_coverage"} if flagged else set()
        metrics = [
            ("ENST1", True, flags, 150, 250),
            ("ENST2", True, set(), 50, 120),
            ("NM_1", True, {"low_exome_coverage"}, 280, 320),
            ("ENST3", False, {"low_exome_coverage"}, 90, 110),
            ("ENST4", True, {"low_exome_coverage"}, None, None),
        ]
        metrics_ht = hl.Table.parallelize(
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
