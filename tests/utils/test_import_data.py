"""Tests for the import functions in ``gnomad_constraint.utils.import_data``."""

import gzip
from pathlib import Path
from typing import List

import hail as hl
import pytest

from gnomad_constraint.resources.constants import (
    GENOMIC_FEATURE_COLUMNS,
    GENOMIC_FEATURE_MISSING,
)
from gnomad_constraint.utils.import_data import (
    check_unique_windows,
    format_dnm_ht,
    format_genomic_features,
    get_dnm_import_args,
    get_genomic_feature_import_args,
    union_dnm_hts,
)


def _feature_ht(element_ids: List[str]) -> hl.Table:
    """
    Return a formatted feature Table with one window per element ID.

    :param element_ids: BED windows, as "chr-start-end".
    :return: Table from ``format_genomic_features``.
    """
    ht = hl.Table.parallelize(
        [
            {"element_id": e, **dict.fromkeys(GENOMIC_FEATURE_COLUMNS, 1.0)}
            for e in element_ids
        ],
        hl.tstruct(
            element_id=hl.tstr, **dict.fromkeys(GENOMIC_FEATURE_COLUMNS, hl.tfloat64)
        ),
    )
    return format_genomic_features(ht)


class TestGetDnmImportArgs:
    """Test the get_dnm_import_args function."""

    def test_reads_published_format(self, tmp_path: Path) -> None:
        """
        Test that the arguments read a gzip compressed CSV with quoted commas.

        The published DNM files are gzip, not bgzip, so Hail needs ``force``, and
        some of their fields have commas inside quotes.
        """
        path = tmp_path / "dnms.csv.gz"
        with gzip.open(path, "wt") as f:
            f.write(
                "Variant,Study,Role\n"
                'chr1:100:A:G,"ASD, 2026",proband\n'
                'chr2:5:C:T,DDD,"proband, affected"\n'
            )
        ht = hl.import_table(str(path), **get_dnm_import_args())
        assert ht.collect() == [
            hl.Struct(Variant="chr1:100:A:G", Study="ASD, 2026", Role="proband"),
            hl.Struct(Variant="chr2:5:C:T", Study="DDD", Role="proband, affected"),
        ]


class TestUnionDnmHts:
    """Test the union_dnm_hts function."""

    def test_union(self) -> None:
        """
        Test that every file's DNMs are kept, with only their ``Variant`` field.

        The files have different extra columns, so each is reduced to ``Variant``
        before the union. A variant in both files is two rows until
        ``format_dnm_ht``.
        """
        asd_ht = hl.Table.parallelize(
            [{"Variant": "chr1:100:A:G", "Role": "proband"}],
            hl.tstruct(Variant=hl.tstr, Role=hl.tstr),
        )
        ddd_ht = hl.Table.parallelize(
            [{"Variant": v, "LCR": False} for v in ["chr1:100:A:G", "chr2:5:C:T"]],
            hl.tstruct(Variant=hl.tstr, LCR=hl.tbool),
        )
        ht = union_dnm_hts([asd_ht, ddd_ht])
        assert ht.row.dtype == hl.tstruct(Variant=hl.tstr)
        assert sorted(ht.Variant.collect()) == [
            "chr1:100:A:G",
            "chr1:100:A:G",
            "chr2:5:C:T",
        ]


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


class TestGetGenomicFeatureImportArgs:
    """Test the get_genomic_feature_import_args function."""

    def test_spellings_are_missing(self, tmp_path: Path) -> None:
        """
        Test that each missing-value spelling is imported as missing, not NaN.

        Regression test: this test called ``hl.import_table`` with its own
        arguments, so it still passed with ``missing`` dropped from the import in
        ``prepare_data.py --import-genomic-features``.
        """
        n = len(GENOMIC_FEATURE_COLUMNS)
        spellings = (GENOMIC_FEATURE_MISSING * n)[:n]
        path = tmp_path / "features.tsv"
        path.write_text(
            "\t".join(["element_id", *GENOMIC_FEATURE_COLUMNS])
            + "\n"
            + "\t".join(["chr1-0-1000", *spellings])
            + "\n"
            + "\t".join(["chr1-1000-2000", *[str(i / 2) for i in range(n)]])
            + "\n"
        )
        ht = hl.import_table(str(path), **get_genomic_feature_import_args())
        rows = {
            r.element_id: [r[c] for c in GENOMIC_FEATURE_COLUMNS] for r in ht.collect()
        }
        assert rows == {
            "chr1-0-1000": [None] * n,
            "chr1-1000-2000": [i / 2 for i in range(n)],
        }


class TestCheckUniqueWindows:
    """Test the check_unique_windows function."""

    def test_unique(self) -> None:
        """Test that a matrix with each window once passes."""
        check_unique_windows(_feature_ht(["chr1-0-1000", "chr1-1000-2000"]))

    def test_duplicated(self) -> None:
        """Test that a window listed twice raises, with the number of extra rows."""
        ht = _feature_ht(["chr1-0-1000", "chr1-1000-2000", "chr1-0-1000"])
        with pytest.raises(ValueError, match="^1 windows are duplicated"):
            check_unique_windows(ht)
