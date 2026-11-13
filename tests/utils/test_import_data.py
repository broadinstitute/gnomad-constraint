"""Tests for the import functions in ``gnomad_constraint.utils.import_data``."""

from pathlib import Path

import hail as hl

from gnomad_constraint.resources.constants import (
    GENOMIC_FEATURE_COLUMNS,
    GENOMIC_FEATURE_MISSING,
)
from gnomad_constraint.utils.import_data import format_dnm_ht, format_genomic_features


class TestFormatDnmHt:
    """Test the format_dnm_ht function."""

    def test_variants(self) -> None:
        """
        Test that DNMs are keyed by GRCh38 variant, with one row per variant.

        A variant reported in both files is one row. Indels are kept: they're
        removed when the DNMs are joined to the context Table.
        """
        ht = hl.Table.parallelize(
            [{"Variant": v} for v in ["chr1:100:A:G", "chr1:200:C:CT", "chr1:100:A:G"]],
            hl.tstruct(Variant=hl.tstr),
        )
        ht = format_dnm_ht(ht)
        assert ht.row.dtype == hl.tstruct(
            locus=hl.tlocus("GRCh38"), alleles=hl.tarray(hl.tstr)
        )
        assert list(ht.key) == ["locus", "alleles"]
        assert ht.collect() == [
            hl.Struct(locus=hl.Locus("chr1", 100, "GRCh38"), alleles=["A", "G"]),
            hl.Struct(locus=hl.Locus("chr1", 200, "GRCh38"), alleles=["C", "CT"]),
        ]


class TestFormatGenomicFeatures:
    """Test the format_genomic_features function."""

    def test_windows(self) -> None:
        """
        Test that each window is its BED window in 1-based positions.

        The BED window ``chr1-0-1000`` covers positions 1-1000 and ``chr1-1000-2000``
        positions 1001-2000.
        """
        ht = hl.Table.parallelize(
            [
                {
                    "element_id": "chr1-0-1000",
                    **dict.fromkeys(GENOMIC_FEATURE_COLUMNS, 1.0),
                },
                {
                    "element_id": "chr1-1000-2000",
                    **dict.fromkeys(GENOMIC_FEATURE_COLUMNS, 2.0),
                },
            ],
            hl.tstruct(
                element_id=hl.tstr,
                **dict.fromkeys(GENOMIC_FEATURE_COLUMNS, hl.tfloat64),
            ),
        )
        ht = format_genomic_features(ht)
        assert list(ht.row) == ["interval", "element_id", *GENOMIC_FEATURE_COLUMNS]
        loci_ht = hl.Table.parallelize(
            [{"pos": p} for p in [1, 1000, 1001, 2000, 2001]], hl.tstruct(pos=hl.tint32)
        )
        window = ht[hl.locus("chr1", loci_ht.pos, "GRCh38")]
        loci_ht = loci_ht.annotate(
            element_id=window.element_id, dist2telo_1k=window.dist2telo_1k
        )
        rows = {r.pos: (r.element_id, r.dist2telo_1k) for r in loci_ht.collect()}
        assert rows == {
            1: ("chr1-0-1000", 1.0),
            1000: ("chr1-0-1000", 1.0),
            1001: ("chr1-1000-2000", 2.0),
            2000: ("chr1-1000-2000", 2.0),
            2001: (None, None),
        }


class TestGenomicFeatureMissing:
    """Test the GENOMIC_FEATURE_MISSING constant used to import the feature matrix."""

    def test_spellings_are_missing(self, tmp_path: Path) -> None:
        """
        Test that each missing-value spelling is imported as missing, not NaN.

        ``prepare_data.py --import-genomic-features`` imports the matrix with float
        types and ``missing=GENOMIC_FEATURE_MISSING``, so NaN never reaches the rr
        models.
        """
        columns = GENOMIC_FEATURE_COLUMNS[:4]
        path = tmp_path / "features.tsv"
        path.write_text(
            "\t".join(["element_id", *columns])
            + "\n"
            + "\t".join(["chr1-0-1000", "NA", "nan", "NaN", ""])
            + "\n"
            + "\t".join(["chr1-1000-2000", "0.5", "1", "2.5", "3"])
            + "\n"
        )
        ht = hl.import_table(
            str(path),
            types={c: hl.tfloat64 for c in columns},
            missing=GENOMIC_FEATURE_MISSING,
        )
        rows = {r.element_id: [r[c] for c in columns] for r in ht.collect()}
        assert rows == {
            "chr1-0-1000": [None, None, None, None],
            "chr1-1000-2000": [0.5, 1.0, 2.5, 3.0],
        }
