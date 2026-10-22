"""Tests for the preprocessing functions in ``gnomad_constraint.utils.constraint``."""

from typing import Any, Dict, List, Optional

import hail as hl
import pytest
from hail.utils.misc import new_temp_file

from gnomad_constraint.utils.constraint import (
    calculate_mu_by_downsampling,
    create_training_set,
    filter_freq_for_constraint,
    get_annotations_for_computing_mu,
    get_build_calibration_model_annotation,
    get_exome_coverage_expr,
    get_exomes_observed_and_possible,
    prepare_context_ht,
    prepare_ht_for_constraint_calculations,
)

FREQ_TYPE = hl.tarray(hl.tstruct(AC=hl.tint32, AF=hl.tfloat64))

GENOMES_FREQ_META = [
    {"group": "adj"},
    {"group": "adj", "pop": "global", "downsampling": "1000"},
    {"group": "adj", "pop": "afr", "downsampling": "1000"},
    {"group": "raw"},
]
"""Genomes frequency metadata, which labels genetic ancestry groups with ``pop``."""

AN_STRATA_META = [
    {"group": "adj"},
    {"group": "adj", "sex": "XX"},
    {"group": "adj", "sex": "XY"},
]
AN_STRATA_SAMPLE_COUNT = [100, 60, 40]
"""Exome AN strata: 100 samples, 60 XX and 40 XY."""


def _locus(contig: str, pos: int) -> hl.Locus:
    """Return a GRCh38 locus."""
    return hl.Locus(contig, pos, reference_genome="GRCh38")


def _freq(ac: int, af: float, n: int = 2) -> List[hl.Struct]:
    """Return a frequency array with ``n`` identical entries."""
    return [hl.Struct(AC=ac, AF=af)] * n


class TestFilterFreqForConstraint:
    """Test the filter_freq_for_constraint function."""

    FREQ_META = [
        {"group": "adj"},
        {"group": "raw"},
        {"group": "adj", "gen_anc": "afr"},
        {"group": "adj", "gen_anc": "nfe"},
        {"group": "adj", "gen_anc": "global", "downsampling": "100"},
        {"group": "adj", "gen_anc": "afr", "downsampling": "100"},
    ]

    @pytest.mark.parametrize(
        "kwargs, expected_idx",
        [
            ({}, [0]),
            ({"gen_ancs": ["afr"]}, [0, 2]),
            ({"downsamplings": [100]}, [0, 4]),
            (
                {
                    "gen_ancs": ["afr"],
                    "downsamplings": [100],
                    "downsampling_gen_ancs": ["afr"],
                },
                [0, 2, 4, 5],
            ),
            ({"gen_ancs": ["amr"]}, [0]),
        ],
    )
    def test_keeps_requested_entries_in_order(
        self, kwargs: Dict[str, Any], expected_idx: List[int]
    ) -> None:
        """Test that adj is always kept, followed by the requested entries in order."""
        freq_expr = hl.literal(
            [hl.Struct(AC=i, AF=0.0) for i in range(len(self.FREQ_META))]
        )
        freq_expr, meta = filter_freq_for_constraint(
            freq_expr, self.FREQ_META, **kwargs
        )
        assert meta == [self.FREQ_META[i] for i in expected_idx]
        assert [x.AC for x in hl.eval(freq_expr)] == expected_idx


class TestGetAnnotationsForComputingMu:
    """Test the get_annotations_for_computing_mu function."""

    # name, locus, filters, (adj AC, global downsampling AC), mean coverage, GERP,
    # most severe consequence.
    ROWS = [
        ("kept", ("chr1", 1000), set(), (3, 1), 30, 0.0, "intron_variant"),
        ("ds_unobserved", ("chr1", 1001), set(), (3, 0), 30, 0.0, "intergenic_variant"),
        ("absent", ("chr1", 1002), None, None, 30, 0.0, "intron_variant"),
        ("chrx_nonpar", ("chrX", 10000000), set(), (3, 1), 30, 0.0, "intron_variant"),
        ("chrx_par", ("chrX", 100000), set(), (3, 1), 30, 0.0, "intron_variant"),
        ("cov_14", ("chr1", 1003), set(), (3, 1), 14, 0.0, "intron_variant"),
        ("cov_15", ("chr1", 1004), set(), (3, 1), 15, 0.0, "intron_variant"),
        ("cov_60", ("chr1", 1005), set(), (3, 1), 60, 0.0, "intron_variant"),
        ("cov_61", ("chr1", 1006), set(), (3, 1), 61, 0.0, "intron_variant"),
        ("cov_missing", ("chr1", 1007), set(), (3, 1), None, 0.0, "intron_variant"),
        ("gerp_at_lower", ("chr1", 1008), set(), (3, 1), 30, -3.9885, "intron_variant"),
        ("gerp_at_upper", ("chr1", 1009), set(), (3, 1), 30, 2.6607, "intron_variant"),
        ("gerp_inside", ("chr1", 1010), set(), (3, 1), 30, 2.66, "intron_variant"),
        ("missense", ("chr1", 1011), set(), (3, 1), 30, 0.0, "missense_variant"),
        ("filtered", ("chr1", 1012), {"AC0"}, (3, 1), 30, 0.0, "intron_variant"),
        ("ds_ac_5", ("chr1", 1013), set(), (9, 5), 30, 0.0, "intron_variant"),
        ("ds_ac_6", ("chr1", 1014), set(), (9, 6), 30, 0.0, "intron_variant"),
    ]

    def _ht(self, gen_ancs: Optional[List[str]] = None) -> hl.Table:
        """Return the sites in ``ROWS`` annotated with the mutation rate annotations."""
        rows = []
        for name, (contig, pos), filters, ac, cov, gerp, csq in self.ROWS:
            # Entries follow GENOMES_FREQ_META: adj, global ds, afr ds, raw.
            freq = (
                None
                if ac is None
                else [
                    hl.Struct(AC=ac[0], AF=ac[0] / 1e5),
                    hl.Struct(AC=ac[1], AF=ac[1] / 2e3),
                    hl.Struct(AC=ac[1], AF=ac[1] / 2e3),
                    hl.Struct(AC=ac[0], AF=ac[0] / 1e5),
                ]
            )
            rows.append(
                dict(
                    name=name,
                    locus=_locus(contig, pos),
                    filters=filters,
                    freq=freq,
                    cov=cov,
                    gerp=gerp,
                    csq=csq,
                )
            )
        ht = hl.Table.parallelize(
            rows,
            hl.tstruct(
                name=hl.tstr,
                locus=hl.tlocus("GRCh38"),
                filters=hl.tset(hl.tstr),
                freq=FREQ_TYPE,
                cov=hl.tint32,
                gerp=hl.tfloat64,
                csq=hl.tstr,
            ),
            key="name",
        )
        expr, globals_expr = get_annotations_for_computing_mu(
            ht.locus,
            ht.filters,
            ht.freq,
            GENOMES_FREQ_META,
            ht.cov,
            ht.gerp,
            ht.csq,
            gen_ancs=gen_ancs,
        )
        return ht.annotate(**expr).annotate_globals(calculate_mu_globals=globals_expr)

    @pytest.fixture(scope="class")
    def mu_rows(self) -> Dict[str, hl.Struct]:
        """Return the annotated sites by name, collected once for the class."""
        return {r.name: r for r in self._ht().collect()}

    @pytest.mark.parametrize(
        "name, observed, possible",
        [
            ("kept", [1, 1], 1),
            ("ds_unobserved", [1, 0], 1),
            ("chrx_nonpar", None, None),
            ("chrx_par", None, None),
            ("cov_14", None, None),
            ("cov_15", [1, 1], 1),
            ("cov_60", [1, 1], 1),
            ("cov_61", None, None),
            ("cov_missing", None, None),
            ("gerp_at_lower", None, None),
            ("gerp_at_upper", None, None),
            ("gerp_inside", [1, 1], 1),
            ("missense", None, None),
            ("filtered", None, None),
            ("ds_ac_5", [1, 1], 1),
            ("ds_ac_6", None, None),
        ],
    )
    def test_observed_and_possible(
        self,
        mu_rows: Dict[str, hl.Struct],
        name: str,
        observed: Optional[List[int]],
        possible: Optional[int],
    ) -> None:
        """Test which sites are kept for the mutation rate, and their counts.

        Sites must be autosomal (PAR excluded), have mean genome coverage in
        [15, 60], GERP strictly between the cutoffs, an intron or intergenic most
        severe consequence, and pass genome filters with AC <= 5 at the
        downsampling. Excluded sites get missing counts.
        """
        row = mu_rows[name]
        assert row.observed_variants == observed
        assert row.possible_variants == possible

    def test_site_absent_from_genomes(self, mu_rows: Dict[str, hl.Struct]) -> None:
        """Test that a site absent from the genomes is possible but not observed.

        The filtered frequency array of an absent site holds missing entries rather
        than being missing itself, so each observed count is 0.
        """
        row = mu_rows["absent"]
        assert row.genomes_freq == [None, None]
        assert row.observed_variants == [0, 0]
        assert row.possible_variants == 1

    def test_genomes_freq_is_filtered(self, mu_rows: Dict[str, hl.Struct]) -> None:
        """Test that the returned frequency array keeps adj and the downsampling."""
        row = mu_rows["ds_ac_5"]
        assert [f.AC for f in row.genomes_freq] == [9, 5]

    @pytest.mark.parametrize(
        "gen_ancs, freq_meta",
        [
            (None, GENOMES_FREQ_META[:2]),
            (["afr"], GENOMES_FREQ_META[:3]),
        ],
    )
    def test_globals(
        self, gen_ancs: Optional[List[str]], freq_meta: List[Dict[str, str]]
    ) -> None:
        """Test the globals, which always include the global downsampling."""
        ht = self._ht(gen_ancs)
        g = hl.eval(ht.calculate_mu_globals)
        assert g.freq_meta == freq_meta
        assert g.downsampling_idx == 1
        assert g.downsampling_level == 1000
        assert g.genetic_ancestry_groups == gen_ancs
        assert (g.ac_cutoff, g.min_cov, g.max_cov) == (5, 15, 60)
        assert (g.gerp_lower_cutoff, g.gerp_upper_cutoff) == (-3.9885, 2.6607)
        assert g.most_severe_consequence == ["intron_variant", "intergenic_variant"]
        row = {r.name: r for r in ht.collect()}["kept"]
        assert len(row.observed_variants) == len(freq_meta)


class TestGetExomeCoverageExpr:
    """Test the get_exome_coverage_expr function."""

    def _ht(self, coverage_fields: Dict[str, Any]) -> hl.Table:
        """Return sites on an autosome, chrX/chrY non-PAR, and the chrX PAR."""
        rows = [
            dict(
                name=name,
                locus=_locus(contig, pos),
                coverage=hl.Struct(exomes=hl.Struct(**coverage_fields)),
                AN=hl.Struct(exomes=an),
            )
            for name, contig, pos, an in [
                ("autosome", "chr1", 1000, 150),
                ("autosome_boundary", "chr1", 1001, 199),
                ("chrx_nonpar", "chrX", 10000000, 80),
                ("chry_nonpar", "chrY", 3000000, 20),
                ("chrx_par", "chrX", 100000, 100),
                ("missing_an", "chr1", 1002, None),
            ]
        ]
        ht = hl.Table.parallelize(
            rows,
            hl.tstruct(
                name=hl.tstr,
                locus=hl.tlocus("GRCh38"),
                coverage=hl.tstruct(
                    exomes=hl.tstruct(**{k: hl.tint32 for k in coverage_fields})
                ),
                AN=hl.tstruct(exomes=hl.tint32),
            ),
            key="name",
        )
        return ht.annotate_globals(
            an_globals=hl.struct(
                exomes=hl.struct(
                    strata_meta=AN_STRATA_META,
                    strata_sample_count=AN_STRATA_SAMPLE_COUNT,
                )
            )
        )

    def _coverage(self, ht: hl.Table, metric: str) -> Dict[str, Optional[int]]:
        """Return the coverage expression's value for each site, by site name."""
        ht = ht.annotate(cov=get_exome_coverage_expr(ht, metric))
        return {r.name: r.cov for r in ht.collect()}

    def test_an_percent(self) -> None:
        """Test AN as a percent of the alleles possible on each chromosome.

        The denominator is 2 * all samples on autosomes and in the PAR,
        2 * XX + XY samples on chrX non-PAR, and XY samples on chrY non-PAR. The
        result is truncated to an integer.
        """
        cov = self._coverage(self._ht({"median_approx": 30}), "AN_percent")
        assert cov == {
            "autosome": 75,
            "autosome_boundary": 99,
            "chrx_nonpar": 50,
            "chry_nonpar": 50,
            "chrx_par": 50,
            "missing_an": None,
        }

    def test_an(self) -> None:
        """Test that the AN metric is the exome allele number."""
        cov = self._coverage(self._ht({"median_approx": 30}), "AN")
        assert cov["autosome"] == 150
        assert cov["missing_an"] is None

    @pytest.mark.parametrize(
        "coverage_fields, expected",
        [
            ({"median_approx": 30, "median": 20}, 30),
            ({"median": 20}, 20),
        ],
    )
    def test_median(self, coverage_fields: Dict[str, int], expected: int) -> None:
        """Test that the median metric prefers median_approx over median."""
        cov = self._coverage(self._ht(coverage_fields), "median")
        assert cov["autosome"] == expected

    def test_invalid_metric(self) -> None:
        """Test that an unknown coverage metric raises an error."""
        with pytest.raises(ValueError, match="Exome coverage metric must be one of"):
            get_exome_coverage_expr(self._ht({"median": 20}), "mean")


class TestGetExomesObservedAndPossible:
    """Test the get_exomes_observed_and_possible function."""

    FREQ_META = [{"group": "adj"}, {"group": "adj", "gen_anc": "afr"}]

    def _ht(self) -> hl.Table:
        """Return sites covering each observed/possible case, with obs/pos annotated."""
        rows = [
            # name, group, filters, freq, coverage
            ("observed", "g1", set(), _freq(1, 1e-5), 30),
            ("absent", "g1", None, None, 30),
            ("filtered", "g2", {"AC0"}, _freq(0, 0.0), 30),
            ("absent2", "g2", None, None, 30),
            ("common", "g3", set(), _freq(500, 0.2), 30),
            ("no_coverage", "g4", None, None, None),
        ]
        ht = hl.Table.parallelize(
            [
                dict(name=n, grp=g, filters=f, freq=fr, coverage=c)
                for n, g, f, fr, c in rows
            ],
            hl.tstruct(
                name=hl.tstr,
                grp=hl.tstr,
                filters=hl.tset(hl.tstr),
                freq=FREQ_TYPE,
                coverage=hl.tint32,
            ),
            key="name",
        )
        expr, _ = get_exomes_observed_and_possible(
            ht.filters, ht.freq, self.FREQ_META, ht.coverage, gen_ancs=["afr"]
        )
        return ht.annotate(**expr)

    @pytest.mark.parametrize(
        "name, observed, possible",
        [
            ("observed", [1, 1], 1),
            ("absent", [0, 0], 1),
            ("filtered", [0, 0], 1),
            ("common", [0, 0], 0),
            ("no_coverage", None, None),
        ],
    )
    def test_per_variant_counts(
        self, name: str, observed: Optional[List[int]], possible: Optional[int]
    ) -> None:
        """Test observed and possible counts for each kind of site."""
        row = {r.name: r for r in self._ht().collect()}[name]
        assert row.observed_variants == observed
        assert row.possible_variants == possible

    def test_unobserved_group_sums_to_zero(self) -> None:
        """Test that a group with no observed variants sums to 0, not missing.

        Regression test for gnomad_production#2294: unobserved possible variants
        used to get a missing observed count, so a group with no observed
        variants summed to missing.
        """
        ht = self._ht()
        ht = ht.group_by("grp").aggregate(
            obs=hl.agg.array_sum(ht.observed_variants),
            pos=hl.agg.sum(ht.possible_variants),
        )
        sums = {r.grp: (r.obs, r.pos) for r in ht.collect()}
        assert sums["g1"] == ([1, 1], 2)
        assert sums["g2"] == ([0, 0], 2)
        assert sums["g3"] == ([0, 0], 0)


class TestGetBuildCalibrationModelAnnotation:
    """Test the get_build_calibration_model_annotation function."""

    CSQ_TYPE = hl.tstruct(
        transcript_id=hl.tstr,
        most_severe_consequence=hl.tstr,
        canonical=hl.tint32,
        mane_select=hl.tstr,
    )

    # name, (transcript ID, most severe consequence, canonical, MANE Select),
    # exome coverage, CpG.
    ROWS = [
        ("mane_syn_high", ("ENST1", "synonymous_variant", 1, "NM_1"), 50, False),
        ("mane_syn_low_cpg", ("ENST1", "synonymous_variant", 1, "NM_1"), 20, True),
        ("canonical_syn", ("ENST2", "synonymous_variant", 1, None), 50, False),
        ("mane_missense", ("ENST1", "missense_variant", 1, "NM_1"), 50, False),
        ("refseq_syn", ("NM_1", "synonymous_variant", 1, "NM_1"), 50, False),
        ("mane_syn_zero_cov", ("ENST1", "synonymous_variant", 1, "NM_1"), 0, False),
        ("mane_syn_cov_70", ("ENST1", "synonymous_variant", 1, "NM_1"), 70, False),
        ("no_csqs", None, 50, False),
    ]

    def _build_models(self, **kwargs: Any) -> Dict[str, Optional[hl.Struct]]:
        """Return the build model annotation for each site in ``ROWS``, by name."""
        rows = [
            dict(
                name=name,
                csqs=(
                    []
                    if csq is None
                    else [
                        hl.Struct(
                            transcript_id=csq[0],
                            most_severe_consequence=csq[1],
                            canonical=csq[2],
                            mane_select=csq[3],
                        )
                    ]
                ),
                cov=cov,
                cpg=cpg,
            )
            for name, csq, cov, cpg in self.ROWS
        ]
        ht = hl.Table.parallelize(
            rows,
            hl.tstruct(
                name=hl.tstr,
                csqs=hl.tarray(self.CSQ_TYPE),
                cov=hl.tint32,
                cpg=hl.tbool,
            ),
            key="name",
        )
        ht = ht.annotate(
            build_model=get_build_calibration_model_annotation(
                ht.cov, ht.csqs, ht.cpg, hl.literal("autosome_or_par"), **kwargs
            )
        )
        return {r.name: r.build_model for r in ht.collect()}

    def test_mane_select(self) -> None:
        """Test model groups when filtering to MANE Select synonymous variants.

        Only Ensembl MANE Select synonymous variants get a group. CpG is kept only
        in the high coverage model, and coverage 0 is in neither model.
        """
        bm = self._build_models()
        assert bm["mane_syn_high"] == hl.Struct(
            high_or_low_coverage="high",
            model_group=hl.Struct(cpg=False, genomic_region="autosome_or_par"),
        )
        assert bm["mane_syn_low_cpg"] == hl.Struct(
            high_or_low_coverage="low",
            model_group=hl.Struct(cpg=None, genomic_region="autosome_or_par"),
        )
        assert bm["mane_syn_cov_70"].high_or_low_coverage == "high"
        for name in [
            "canonical_syn",
            "mane_missense",
            "refseq_syn",
            "mane_syn_zero_cov",
            "no_csqs",
        ]:
            assert bm[name] is None, name

    def test_canonical(self) -> None:
        """Test that canonical mode keeps canonical transcripts without MANE Select."""
        bm = self._build_models(synonymous_transcript_filter_field="canonical")
        assert bm["canonical_syn"].high_or_low_coverage == "high"
        assert bm["mane_syn_high"].high_or_low_coverage == "high"
        assert bm["mane_missense"] is None

    def test_skip_coverage_model(self) -> None:
        """Test that low coverage sites get no group without a coverage model."""
        bm = self._build_models(skip_coverage_model=True)
        assert bm["mane_syn_high"].high_or_low_coverage == "high"
        assert bm["mane_syn_low_cpg"] is None

    def test_coverage_cutoffs(self) -> None:
        """Test the low and upper coverage cutoffs."""
        bm = self._build_models(low_cov_cutoff=20, upper_cov_cutoff=60)
        assert bm["mane_syn_low_cpg"] is None
        assert bm["mane_syn_cov_70"] is None
        assert bm["mane_syn_high"].high_or_low_coverage == "high"

    def test_invalid_filter_field(self) -> None:
        """Test that an unknown transcript filter field raises an error."""
        csqs = hl.empty_array(self.CSQ_TYPE)
        with pytest.raises(ValueError, match="must be either 'canonical' or"):
            get_build_calibration_model_annotation(
                hl.int32(50),
                csqs,
                hl.bool(False),
                hl.literal("autosome_or_par"),
                synonymous_transcript_filter_field="ensembl",
            )


class TestCreateTrainingSet:
    """Test the create_training_set function."""

    BUILD_MODEL_TYPE = hl.tstruct(
        high_or_low_coverage=hl.tstr,
        model_group=hl.tstruct(cpg=hl.tbool, genomic_region=hl.tstr),
    )

    def _ht(self) -> hl.Table:
        """Return preprocessed sites with build model groups."""
        high = hl.Struct(
            high_or_low_coverage="high",
            model_group=hl.Struct(cpg=True, genomic_region="autosome_or_par"),
        )
        low = hl.Struct(
            high_or_low_coverage="low",
            model_group=hl.Struct(cpg=None, genomic_region="autosome_or_par"),
        )
        # name, build model, coverage, observed, possible.
        sites = [
            ("a", high, 50, [1, 0], 1),
            ("b", high, 50, [0, 0], 1),
            ("c", high, 50, [1, 1], 1),
            ("d", low, 20, [1, 1], 1),
            ("no_build_model", None, 50, [1, 1], 1),
            ("not_possible", high, 50, [0, 0], 0),
        ]
        rows = [
            dict(
                name=name,
                context="ACG",
                ref="C",
                alt="T",
                methylation_level=2,
                cpg=True,
                transition=True,
                mutation_type="CpG",
                mutation_type_model="CpG",
                genomic_region="autosome_or_par",
                exomes_coverage=cov,
                calibrate_mu=hl.Struct(
                    observed_variants=obs, possible_variants=pos, build_model=bm
                ),
            )
            for name, bm, cov, obs, pos in sites
        ]
        return hl.Table.parallelize(
            rows,
            hl.tstruct(
                name=hl.tstr,
                context=hl.tstr,
                ref=hl.tstr,
                alt=hl.tstr,
                methylation_level=hl.tint32,
                cpg=hl.tbool,
                transition=hl.tbool,
                mutation_type=hl.tstr,
                mutation_type_model=hl.tstr,
                genomic_region=hl.tstr,
                exomes_coverage=hl.tint32,
                calibrate_mu=hl.tstruct(
                    observed_variants=hl.tarray(hl.tint32),
                    possible_variants=hl.tint32,
                    build_model=self.BUILD_MODEL_TYPE,
                ),
            ),
            key="name",
        )

    def _mutation_ht(self) -> hl.Table:
        """Return a one-row mutation rate Table for ACG C>T at methylation level 2."""
        return hl.Table.parallelize(
            [dict(context="ACG", ref="C", alt="T", methylation_level=2, mu_snp=1e-7)],
            hl.tstruct(
                context=hl.tstr,
                ref=hl.tstr,
                alt=hl.tstr,
                methylation_level=hl.tint32,
                mu_snp=hl.tfloat64,
            ),
            key=["context", "ref", "alt", "methylation_level"],
        )

    def test_counts_by_group(self) -> None:
        """Test that counts are summed per model group and coverage.

        Sites without a build model group or with no possible variants are
        dropped, and every group gets its mutation rate.
        """
        ht = create_training_set(self._ht(), self._mutation_ht(), partition_hint=1)
        rows = {
            (r.build_model.high_or_low_coverage, r.exomes_coverage): r
            for r in ht.collect()
        }
        assert set(rows) == {("high", 50), ("low", 20)}
        assert rows[("high", 50)].observed_variants == [2, 1]
        assert rows[("high", 50)].possible_variants == 3
        assert rows[("low", 20)].observed_variants == [1, 1]
        assert rows[("low", 20)].possible_variants == 1
        assert all(r.mu_snp == pytest.approx(1e-7) for r in rows.values())


class TestCalculateMuByDownsampling:
    """Test the calculate_mu_by_downsampling function."""

    def _ht(self) -> hl.Table:
        """Return sites in two substitution groups with per-group observed counts."""
        # context, ref, alt, observed (adj, downsampling), possible.
        sites = [
            ("ACG", "C", "T", [1, 1], 1),
            ("ACG", "C", "T", [1, 0], 1),
            ("ACG", "C", "T", [0, 0], 1),
            ("ACG", "C", "T", [0, 0], 1),
            ("ACA", "C", "A", [1, 0], 1),
            ("ACA", "C", "A", [0, 0], 1),
            ("ACA", "C", "A", None, None),
        ]
        rows = [
            dict(
                idx=i,
                context=c,
                ref=r,
                alt=a,
                methylation_level=0,
                compute_mu=hl.Struct(observed_variants=obs, possible_variants=pos),
            )
            for i, (c, r, a, obs, pos) in enumerate(sites)
        ]
        ht = hl.Table.parallelize(
            rows,
            hl.tstruct(
                idx=hl.tint32,
                context=hl.tstr,
                ref=hl.tstr,
                alt=hl.tstr,
                methylation_level=hl.tint32,
                compute_mu=hl.tstruct(
                    observed_variants=hl.tarray(hl.tint32),
                    possible_variants=hl.tint32,
                ),
            ),
            key="idx",
        )
        return ht.annotate_globals(calculate_mu_globals=hl.struct(downsampling_idx=1))

    def test_mutation_rate(self) -> None:
        """Test the proportion observed and mutation rate formula.

        There are 6 possible variants, so 2 bases. Summed observed counts are
        [3, 1], giving correction factors of 1.2e-8 / ([3, 1] / 2). ``mu_snp`` is
        the mutation rate at the downsampling index.
        """
        ht = calculate_mu_by_downsampling(self._ht())
        rows = {(r.context, r.ref, r.alt): r for r in ht.collect()}
        cpg, tv = rows[("ACG", "C", "T")], rows[("ACA", "C", "A")]
        correction = [1.2e-8 / (3 / 2), 1.2e-8 / (1 / 2)]

        assert cpg.possible_variants == 4 and cpg.observed_variants == [2, 1]
        assert cpg.proportion_observed == pytest.approx([0.5, 0.25])
        assert cpg.mu == pytest.approx([correction[0] * 0.5, correction[1] * 0.25])
        assert cpg.mu_snp == pytest.approx(correction[1] * 0.25)

        assert tv.possible_variants == 2 and tv.observed_variants == [1, 0]
        assert tv.mu == pytest.approx([correction[0] * 0.5, 0.0])
        assert tv.mu_snp == pytest.approx(0.0)

    def test_mutation_type(self) -> None:
        """Test that mutation type annotations are added."""
        ht = calculate_mu_by_downsampling(self._ht())
        rows = {(r.context, r.ref, r.alt): r for r in ht.collect()}
        assert rows[("ACG", "C", "T")].mutation_type == "CpG"
        assert rows[("ACG", "C", "T")].cpg
        assert rows[("ACA", "C", "A")].mutation_type == "transversion"
        assert not rows[("ACA", "C", "A")].transition


CSQ_TYPE = hl.tstruct(
    transcript_id=hl.tstr,
    gene_id=hl.tstr,
    gene_symbol=hl.tstr,
    biotype=hl.tstr,
    consequence_terms=hl.tarray(hl.tstr),
    mane_select=hl.tstr,
    canonical=hl.tint32,
    lof=hl.tstr,
    lof_filter=hl.tstr,
    lof_flags=hl.tstr,
    lof_info=hl.tstr,
)
"""VEP transcript consequence fields used by prepare_context_ht."""

# name, locus, heptamer context, consequence, LOFTEE (lof, lof_filter, lof_info),
# exomes AF (missing when absent from the exomes), exomes filters.
CONTEXT_SITES = [
    ("syn_cpg", ("chr1", 1000), "AAACGAA", "synonymous_variant", None, 1e-5, set()),
    (
        "lof_end_trunc",
        ("chr1", 2000),
        "AAACAAA",
        "stop_gained",
        ("HC", None, "GERP_DIST:-5.0,50_BP_RULE:FAIL"),
        None,
        None,
    ),
    (
        "lof_hc",
        ("chr1", 3000),
        "AAACAAA",
        "stop_gained",
        ("HC", None, "GERP_DIST:5.0,50_BP_RULE:PASS"),
        0.0,
        {"AC0"},
    ),
    ("chrx_nonpar", ("chrX", 10000000), "AAACAAA", "intron_variant", None, 0.3, set()),
    ("chrx_par", ("chrX", 100000), "AAACAAA", "intron_variant", None, 0.7, set()),
    ("chry_nonpar", ("chrY", 3000000), "AAACAAA", "intron_variant", None, None, None),
]


def _context_inputs() -> Dict[str, Any]:
    """Return the context Table and the resource Tables prepare_context_ht joins."""
    loci = {name: _locus(*locus) for name, locus, *_ in CONTEXT_SITES}

    def _csq(csq: str, lof: Optional[tuple]) -> hl.Struct:
        lof = lof or (None, None, None)
        return hl.Struct(
            transcript_id="ENST1",
            gene_id="ENSG1",
            gene_symbol="G1",
            biotype="protein_coding",
            consequence_terms=[csq],
            mane_select="NM_1",
            canonical=1,
            lof=lof[0],
            lof_filter=lof[1],
            lof_flags=None,
            lof_info=lof[2],
        )

    context_ht = hl.Table.parallelize(
        [
            dict(
                locus=loci[name],
                alleles=["C", "T"],
                context=context,
                was_split=False,
                vep=hl.Struct(
                    most_severe_consequence=csq,
                    transcript_consequences=[_csq(csq, lof)],
                ),
            )
            for name, _, context, csq, lof, *_ in CONTEXT_SITES
        ],
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            context=hl.tstr,
            was_split=hl.tbool,
            vep=hl.tstruct(
                most_severe_consequence=hl.tstr,
                transcript_consequences=hl.tarray(CSQ_TYPE),
            ),
        ),
        key=["locus", "alleles"],
    )

    def _locus_ht(field_values: Dict[str, Any], dtypes: Dict[str, hl.HailType]):
        return hl.Table.parallelize(
            [dict(locus=loci[n], **field_values[n]) for n in field_values],
            hl.tstruct(locus=hl.tlocus("GRCh38"), **dtypes),
            key="locus",
        )

    def _coverage_ht(mean: int) -> hl.Table:
        return _locus_ht(
            {n: dict(mean=mean, median_approx=mean) for n in loci},
            dict(mean=hl.tint32, median_approx=hl.tint32),
        )

    def _an_ht(an: int) -> hl.Table:
        ht = _locus_ht(
            {n: dict(AN=[an, an]) for n in loci}, dict(AN=hl.tarray(hl.tint32))
        )
        return ht.annotate_globals(
            strata_meta=AN_STRATA_META, strata_sample_count=AN_STRATA_SAMPLE_COUNT
        )

    exomes = {name: (af, filters) for name, _, _, _, _, af, filters in CONTEXT_SITES}
    exomes_freq_ht = hl.Table.parallelize(
        [
            dict(
                locus=loci[name],
                alleles=["C", "T"],
                freq=[hl.Struct(AC=int(af * 1e5), AF=af)] * 2,
                filters=filters,
            )
            for name, (af, filters) in exomes.items()
            if af is not None
        ],
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            freq=FREQ_TYPE,
            filters=hl.tset(hl.tstr),
        ),
        key=["locus", "alleles"],
    ).annotate_globals(freq_meta=[{"group": "adj"}, {"group": "raw"}])
    genomes_freq_ht = hl.Table.parallelize(
        [
            dict(
                locus=loci["syn_cpg"],
                alleles=["C", "T"],
                freq=[hl.Struct(AC=2, AF=1e-4)] * 4,
                filters=set(),
            )
        ],
        exomes_freq_ht.row.dtype,
        key=["locus", "alleles"],
    ).annotate_globals(freq_meta=GENOMES_FREQ_META)

    adj_r = {
        n: dict(adj_r={"ACG": 0.9, "ACA": 1.1}) for n in loci if "chr1" in str(loci[n])
    }
    return dict(
        ht=context_ht,
        coverage_hts={"exomes": _coverage_ht(50), "genomes": _coverage_ht(30)},
        an_hts={"exomes": _an_ht(150), "genomes": _an_ht(100)},
        freq_hts={
            "exomes": exomes_freq_ht.select("freq"),
            "genomes": genomes_freq_ht.select("freq"),
        },
        filter_hts={
            "exomes": exomes_freq_ht.select("filters"),
            "genomes": genomes_freq_ht.select("filters"),
        },
        methylation_ht=_locus_ht(
            {n: dict(methylation_level=10) for n in loci},
            dict(methylation_level=hl.tint32),
        ),
        gerp_ht=_locus_ht({n: dict(S=0.5) for n in loci}, dict(S=hl.tfloat64)),
        adj_r_ht=_locus_ht(adj_r, dict(adj_r=hl.tdict(hl.tstr, hl.tfloat64))),
        syn_adj_r_ht=_locus_ht(
            {n: dict(adj_r={"ACG": 0.8, "ACA": 1.0}) for n in adj_r},
            dict(adj_r=hl.tdict(hl.tstr, hl.tfloat64)),
        ),
    )


@pytest.fixture(scope="module")
def annotated_context_ht() -> hl.Table:
    """Return the output of prepare_context_ht, checkpointed once for the module."""
    ht = prepare_context_ht(**_context_inputs())
    return ht.checkpoint(new_temp_file("prepare_context_ht", "ht"))


def _rows_by_site(ht: hl.Table) -> Dict[str, hl.Struct]:
    """Return the rows of a Table keyed by locus, by ``CONTEXT_SITES`` name."""
    names = {_locus(*locus): name for name, locus, *_ in CONTEXT_SITES}
    return {names[r.locus]: r for r in ht.collect()}


class TestPrepareContextHt:
    """Test the prepare_context_ht function."""

    @pytest.fixture(scope="class")
    def context_rows(self, annotated_context_ht: hl.Table) -> Dict[str, hl.Struct]:
        """Return the annotated context rows by site name."""
        return _rows_by_site(annotated_context_ht)

    def test_genomic_region(self, context_rows: Dict[str, hl.Struct]) -> None:
        """Test the genomic region labels, with the chrX PAR as autosome_or_par."""
        rows = context_rows
        assert {n: r.genomic_region for n, r in rows.items()} == {
            "syn_cpg": "autosome_or_par",
            "lof_end_trunc": "autosome_or_par",
            "lof_hc": "autosome_or_par",
            "chrx_nonpar": "chrx_nonpar",
            "chrx_par": "autosome_or_par",
            "chry_nonpar": "chry_nonpar",
        }

    def test_sfs_bin(self, context_rows: Dict[str, hl.Struct]) -> None:
        """Test the SFS bin from the exomes AF.

        A site absent from the exomes and an AF of 0 are both bin 0, and an AF
        above the last cutoff (0.5) has no bin.
        """
        rows = context_rows
        assert {n: r.sfs_bin for n, r in rows.items()} == {
            "syn_cpg": 4,
            "lof_end_trunc": 0,
            "lof_hc": 0,
            "chrx_nonpar": 8,
            "chrx_par": None,
            "chry_nonpar": 0,
        }

    def test_adj_r(self, context_rows: Dict[str, hl.Struct]) -> None:
        """Test that adj_r and syn_adj_r are looked up by locus and trimer context."""
        rows = context_rows
        assert rows["syn_cpg"].context == "ACG"
        assert rows["syn_cpg"].adj_r == pytest.approx(0.9)
        assert rows["syn_cpg"].syn_adj_r == pytest.approx(0.8)
        assert rows["lof_hc"].adj_r == pytest.approx(1.1)
        assert rows["chrx_nonpar"].adj_r is None

    def test_loftee_end_trunc(self, context_rows: Dict[str, hl.Struct]) -> None:
        """Test that HC LoF at an unconstrained end truncation becomes LC."""
        rows = context_rows
        assert rows["lof_end_trunc"].vep.transcript_consequences[0].lof == "LC"
        assert rows["lof_hc"].vep.transcript_consequences[0].lof == "HC"

    def test_coverage_and_an(self, context_rows: Dict[str, hl.Struct]) -> None:
        """Test that coverage and AN are reshaped per data type."""
        row = context_rows["syn_cpg"]
        assert row.coverage.exomes == hl.Struct(mean=50, median_approx=50)
        assert row.coverage.genomes.mean == 30
        assert row.AN == hl.Struct(exomes=150, genomes=100)


class TestPrepareHtForConstraintCalculations:
    """Test the prepare_ht_for_constraint_calculations function."""

    @pytest.fixture(scope="class")
    def preprocessed(self, annotated_context_ht: hl.Table):
        """Return the preprocessed rows by site name and the evaluated globals.

        The coverage model is skipped so the globals test can check that the
        setting is recorded; it does not change the high coverage site's groups.
        """
        ht = prepare_ht_for_constraint_calculations(
            annotated_context_ht,
            exome_coverage_metric="AN_percent",
            skip_coverage_model=True,
        )
        return _rows_by_site(ht), hl.eval(ht.globals)

    def test_annotations(self, preprocessed) -> None:
        """Test the coverage, observed/possible, and model group annotations."""
        rows, _ = preprocessed
        syn = rows["syn_cpg"]
        assert "freq" not in syn
        assert syn.exomes_coverage == 75
        assert syn.calibrate_mu.observed_variants == [1]
        assert syn.calibrate_mu.possible_variants == 1
        assert syn.calibrate_mu.build_model.high_or_low_coverage == "high"
        assert syn.calibrate_mu.apply_model.high_or_low_coverage == "high"
        # Absent from the exomes, and AC0-filtered: possible but not observed.
        for name in ["lof_end_trunc", "lof_hc"]:
            assert rows[name].calibrate_mu.observed_variants == [0], name
            assert rows[name].calibrate_mu.possible_variants == 1, name
            assert rows[name].calibrate_mu.build_model is None, name
        # Genome coverage 30, GERP 0.5, but a synonymous most severe consequence.
        assert syn.compute_mu.possible_variants is None
        # Intron variant, but not autosomal.
        assert rows["chrx_par"].compute_mu.possible_variants is None

    def test_globals(self, preprocessed) -> None:
        """Test that the globals record the parameters used."""
        _, g = preprocessed
        assert g.build_models_globals == hl.Struct(
            synonymous_transcript_filter_field="mane_select",
            low_cov_cutoff=None,
            high_cov_cutoff=40,
            upper_cov_cutoff=None,
            skip_coverage_model=True,
        )
        assert g.apply_models_globals.skip_coverage_model
        assert g.calculate_mu_globals.downsampling_idx == 1
        assert g.exomes_freq_meta == [{"group": "adj"}]
        assert g.max_af == 0.001
