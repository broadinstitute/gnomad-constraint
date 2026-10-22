"""Tests for the resource utility module."""

from typing import Any, Callable, Dict, List

import gnomad.resources.grch37.gnomad as gnomad_grch37
import gnomad.resources.grch37.reference_data as ref_grch37
import gnomad.resources.grch38.gnomad as gnomad_grch38
import gnomad.resources.grch38.reference_data as ref_grch38
import hail as hl
import pytest
from gnomad.resources.resource_utils import ExpressionResource, TableResource

from gnomad_constraint.resources.resource_utils import (
    check_param_scope,
    filter_for_test,
    get_aggregated_expected,
    get_aggregated_per_variant_expected,
    get_annotated_context_ht,
    get_checkpoint_path,
    get_constraint_data,
    get_constraint_group_ht,
    get_constraint_metrics_dataset,
    get_constraint_root,
    get_coverage_ht,
    get_gene_quality_metrics_ht,
    get_lof_threshold_tsv_path,
    get_logging_path,
    get_methylation_ht,
    get_models,
    get_mutation_ht,
    get_per_variant_expected_dataset,
    get_pre_rank_constraint_metrics_dataset,
    get_preprocessed_ht,
    get_release_constraint_ht,
    get_release_constraint_tsv_path,
    get_release_downsampling_tsv_path,
    get_release_mutation_ht,
    get_release_mutation_tsv_path,
    get_sites_resource,
    get_training_dataset,
    get_training_tsv_path,
    get_vep_context_ht,
)

PROD_ROOT = "gs://gnomad/v4.1.1/constraint"
TEST_ROOT = "gs://gnomad-tmp/gnomad_v4.1.1_testing/constraint"


def _path(resource) -> str:
    """Return the path of a resource or of a path string."""
    return resource if isinstance(resource, str) else resource.path


class TestGetConstraintRoot:
    """Test the get_constraint_root function."""

    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({}, "gs://gnomad/v4.1.1/constraint"),
            ({"test": True}, "gs://gnomad-tmp/gnomad_v4.1.1_testing/constraint"),
            ({"temp": True}, "gs://gnomad-tmp/gnomad_v4.1.1/constraint"),
            ({"post_fix": "ukb"}, "gs://gnomad/v4.1.1/constraint_ukb"),
            ({"sub_dir": "metrics"}, "gs://gnomad/v4.1.1/constraint/metrics"),
            (
                {"test": True, "temp": True, "post_fix": "x", "sub_dir": "release"},
                "gs://gnomad-tmp/gnomad_v4.1.1_testing/constraint_x/release",
            ),
        ],
    )
    def test_root(self, kwargs: dict, expected: str) -> None:
        """Test the root for each combination of test, temp, post-fix, and sub_dir."""
        assert get_constraint_root(version="4.1.1", **kwargs) == expected


class TestGetConstraintData:
    """Test the get_constraint_data function."""

    @pytest.mark.parametrize(
        "extension, expected_type",
        [
            ("ht", TableResource),
            ("he", ExpressionResource),
            ("tsv", str),
            ("tsv.bgz", str),
            ("log", str),
        ],
    )
    def test_return_type(self, extension: str, expected_type: type) -> None:
        """Test that the return type depends on the extension."""
        assert isinstance(
            get_constraint_data("x", version="4.1.1", extension=extension),
            expected_type,
        )

    def test_invalid_extension(self) -> None:
        """Test that an unsupported extension raises an error."""
        with pytest.raises(ValueError, match="extension must be one of"):
            get_constraint_data("x", version="4.1.1", extension="csv")

    def test_invalid_version(self) -> None:
        """Test that an unsupported version raises an error."""
        with pytest.raises(ValueError):
            get_constraint_data("x", version="3.1")


class TestReleaseAndGeneQualityPaths:
    """Test the release and gene quality metrics path functions.

    Regression test: these paths used to ignore ``test`` and the post-fixes, so a
    ``--test --overwrite`` run wrote over the production outputs.
    """

    PATHS = [
        (get_release_mutation_ht, "release/gnomad.v4.1.1.mutation_rate{}.ht"),
        (get_release_mutation_tsv_path, "release/gnomad.v4.1.1.mutation_rate{}.tsv"),
        (
            get_release_constraint_ht,
            "release/gnomad.v4.1.1.constraint_metrics{}.ht",
        ),
        (
            get_release_constraint_tsv_path,
            "release/gnomad.v4.1.1.constraint_metrics{}.tsv.bgz",
        ),
        (
            get_release_downsampling_tsv_path,
            "release/gnomad.v4.1.1.constraint_metrics.downsampling{}.tsv.bgz",
        ),
        (
            get_lof_threshold_tsv_path,
            "release/gnomad.v4.1.1.loeuf_percentile_thresholds{}.tsv",
        ),
        (
            get_gene_quality_metrics_ht,
            "metrics/gnomad.v4.1.1.gene_quality_metrics{}.ht",
        ),
    ]

    @pytest.mark.parametrize("func, path", PATHS)
    def test_production_path(self, func: Callable, path: str) -> None:
        """Test the production path."""
        assert _path(func(version="4.1.1")) == (
            "gs://gnomad/v4.1.1/constraint/" + path.format("")
        )

    @pytest.mark.parametrize("func, path", PATHS)
    def test_test_path(self, func: Callable, path: str) -> None:
        """Test that test runs write under the test root."""
        assert _path(func(version="4.1.1", test=True)) == (
            "gs://gnomad-tmp/gnomad_v4.1.1_testing/constraint/" + path.format("")
        )

    @pytest.mark.parametrize("func, path", PATHS)
    def test_post_fixes(self, func: Callable, path: str) -> None:
        """Test that the directory and path post-fixes are applied."""
        out = func(version="4.1.1", directory_post_fix="ukb", path_post_fix="v2")
        assert _path(out) == "gs://gnomad/v4.1.1/constraint_ukb/" + path.format(".v2")


class TestCheckParamScope:
    """Test the check_param_scope function."""

    @pytest.mark.parametrize(
        "kwargs, message",
        [
            ({"data_type": "transcripts"}, "data_type must be one of"),
            ({"model_type": "linear"}, "model_type must be one of"),
            ({"custom_vep_annotation": "worst_csq"}, "custom_vep_annotation must be"),
            ({"extension": "csv"}, "extension must be one of"),
            ({"version": "3.1"}, "The requested version doesn't exist"),
        ],
    )
    def test_invalid_argument(self, kwargs: Dict[str, str], message: str) -> None:
        """Test that each out-of-scope argument raises an error naming it."""
        with pytest.raises(ValueError, match=message):
            check_param_scope(**kwargs)

    @pytest.mark.parametrize(
        "version, build",
        [
            ("2.1.1", "GRCh37"),
            ("4.0", "GRCh38"),
            ("4.1", "GRCh38"),
            ("4.1.1", "GRCh38"),
            ("4.1.2", "GRCh38"),
        ],
    )
    def test_returns_build(self, version: str, build: str) -> None:
        """Test that a valid version returns its genome build."""
        assert check_param_scope(version=version) == build

    def test_no_version(self) -> None:
        """Test that valid arguments without a version return None."""
        assert (
            check_param_scope(
                model_type="plateau",
                custom_vep_annotation="transcript_consequences",
                extension="ht",
                data_type="exomes",
            )
            is None
        )


class TestGetVepContextHt:
    """Test the get_vep_context_ht function."""

    @pytest.mark.parametrize(
        "version, expected",
        [
            ("2.1.1", ref_grch37.vep_context.versions["85"]),
            ("4.0", ref_grch38.vep_context.versions["105"]),
            ("4.1.1", ref_grch38.vep_context.versions["105"]),
        ],
    )
    def test_version_mapping(self, version: str, expected: TableResource) -> None:
        """Test that v2 uses the VEP 85 context and v4 the VEP 105 context."""
        assert get_vep_context_ht(version).path == expected.path

    def test_invalid_version(self) -> None:
        """Test that a version outside v2.1.1 and v4 raises an error."""
        with pytest.raises(ValueError, match="Not a valid gnomAD version"):
            get_vep_context_ht("3.1")


class TestGetSitesResource:
    """Test the get_sites_resource function."""

    @pytest.mark.parametrize(
        "data_type, version, expected",
        [
            (
                "exomes",
                "2.1.1",
                gnomad_grch37.public_release("exomes").versions["2.1.1"],
            ),
            (
                "genomes",
                "2.1.1",
                gnomad_grch37.public_release("genomes").versions["2.1.1"],
            ),
            ("exomes", "4.0", gnomad_grch38.public_release("exomes").versions["4.0"]),
            ("exomes", "4.1", gnomad_grch38.public_release("exomes").versions["4.1"]),
            (
                "exomes",
                "4.1.1",
                gnomad_grch38.public_release("exomes").versions["4.1"],
            ),
            (
                "exomes",
                "4.1.2",
                gnomad_grch38.public_release("exomes").versions["4.1"],
            ),
            (
                "genomes",
                "4.1.1",
                gnomad_grch38.public_release("genomes").versions["3.1.2"],
            ),
            (
                "genomes",
                "4.0",
                gnomad_grch38.public_release("genomes").versions["3.1.2"],
            ),
        ],
    )
    def test_version_mapping(
        self, data_type: str, version: str, expected: TableResource
    ) -> None:
        """Test the sites release used for each pipeline version.

        v4.1.1 and v4.1.2 use the v4.1 exomes, and every v4 version uses the v3.1.2
        genomes because the v4 genomes have no downsamplings.
        """
        assert get_sites_resource(data_type, version).path == expected.path

    def test_invalid_data_type(self) -> None:
        """Test that an unknown data type raises an error."""
        with pytest.raises(ValueError, match="data_type must be one of"):
            get_sites_resource("transcripts", "4.1.1")


class TestGetCoverageHt:
    """Test the get_coverage_ht function."""

    @pytest.mark.parametrize(
        "data_type, version, expected",
        [
            ("exomes", "2.1.1", gnomad_grch37.coverage("exomes")),
            ("genomes", "2.1.1", gnomad_grch37.coverage("genomes")),
            ("exomes", "4.1.1", gnomad_grch38.coverage("exomes")),
            ("genomes", "4.1.1", gnomad_grch38.coverage("genomes")),
        ],
    )
    def test_build_mapping(self, data_type: str, version: str, expected: Any) -> None:
        """Test that the coverage resource comes from the version's genome build."""
        resource = get_coverage_ht(data_type, version)
        assert resource.path == expected.path
        assert list(resource.versions) == list(expected.versions)


class TestGetMethylationHt:
    """Test the get_methylation_ht function."""

    @pytest.mark.parametrize(
        "build, expected",
        [
            ("GRCh37", ref_grch37.methylation_sites),
            ("GRCh38", ref_grch38.methylation_sites),
        ],
    )
    def test_build_mapping(self, build: str, expected: TableResource) -> None:
        """Test that each build returns its methylation resource."""
        assert get_methylation_ht(build) is expected

    def test_invalid_build(self) -> None:
        """Test that an unknown build raises an error."""
        with pytest.raises(ValueError, match="Build must be one of"):
            get_methylation_ht("hg19")


class TestGetModels:
    """Test the get_models function."""

    @pytest.mark.parametrize("model_type", ["plateau", "coverage"])
    def test_path(self, model_type: str) -> None:
        """Test that models are ExpressionResources under the models directory."""
        resource = get_models(model_type, version="4.1.1", path_post_fix="v2")
        assert isinstance(resource, ExpressionResource)
        assert resource.path == (f"{PROD_ROOT}/models/gnomad.v4.1.1.{model_type}.v2.he")

    def test_invalid_model_type(self) -> None:
        """Test that an unknown model type raises an error."""
        with pytest.raises(ValueError, match="model_type must be one of"):
            get_models("linear", version="4.1.1")


class TestPipelineStepPaths:
    """Test the path helpers for the intermediate pipeline outputs."""

    VEP_FUNCS = [
        (get_per_variant_expected_dataset, "apply_models", "per_variant_expected"),
        (
            get_aggregated_per_variant_expected,
            "apply_models",
            "per_variant_expected.aggregated",
        ),
        (get_aggregated_expected, "apply_models", "aggregated_expected"),
        (get_constraint_group_ht, "apply_models", "constraint_group"),
        (get_constraint_metrics_dataset, "metrics", "constraint_metrics"),
        (
            get_pre_rank_constraint_metrics_dataset,
            "metrics",
            "constraint_metrics_pre_rank",
        ),
    ]

    @pytest.mark.parametrize(
        "func, path",
        [
            (
                get_annotated_context_ht,
                "preprocessed_data/gnomad.v4.1.1.annotated_context.ht",
            ),
            (
                get_preprocessed_ht,
                "preprocessed_data/gnomad.v4.1.1.context.preprocessed.ht",
            ),
            (get_mutation_ht, "mutation_rate/gnomad.v4.1.1.mutation_rate.ht"),
            (
                get_training_dataset,
                "training_data/gnomad.v4.1.1.constraint_training.ht",
            ),
            (
                get_training_tsv_path,
                "training_data/gnomad.v4.1.1.constraint_training.tsv.bgz",
            ),
        ],
    )
    def test_path(self, func: Callable, path: str) -> None:
        """Test the production and test paths of each output."""
        assert _path(func(version="4.1.1")) == f"{PROD_ROOT}/{path}"
        assert _path(func(version="4.1.1", test=True)) == f"{TEST_ROOT}/{path}"

    @pytest.mark.parametrize("func, sub_dir, name", VEP_FUNCS)
    @pytest.mark.parametrize(
        "custom_vep_annotation", ["transcript_consequences", "worst_csq_by_gene"]
    )
    def test_custom_vep_annotation_subdir(
        self, func: Callable, sub_dir: str, name: str, custom_vep_annotation: str
    ) -> None:
        """Test that the VEP annotation adds a subdirectory and post-fixes apply."""
        out = func(
            custom_vep_annotation,
            version="4.1.1",
            directory_post_fix="ukb",
            path_post_fix="v2",
        )
        assert out.path == (
            f"gs://gnomad/v4.1.1/constraint_ukb/{sub_dir}/{custom_vep_annotation}/"
            f"gnomad.v4.1.1.{name}.v2.ht"
        )

    @pytest.mark.parametrize("func, sub_dir, name", VEP_FUNCS)
    def test_invalid_custom_vep_annotation(
        self, func: Callable, sub_dir: str, name: str
    ) -> None:
        """Test that an unknown VEP annotation raises an error."""
        with pytest.raises(ValueError, match="custom_vep_annotation must be"):
            func("worst_csq", version="4.1.1")


class TestGetLoggingPath:
    """Test the get_logging_path function."""

    def test_always_under_test_root(self) -> None:
        """Test that log files always go under the test root."""
        assert get_logging_path("constraint", version="4.1.1") == (
            f"{TEST_ROOT}/logging/gnomad.v4.1.1.constraint.log"
        )


class TestGetCheckpointPath:
    """Test the get_checkpoint_path function."""

    def test_always_under_test_root(self) -> None:
        """Test that checkpoints always go under the test root."""
        resource = get_checkpoint_path("tmp", version="4.1.1", directory_post_fix="x")
        assert isinstance(resource, TableResource)
        assert resource.path == (
            "gs://gnomad-tmp/gnomad_v4.1.1_testing/constraint_x/checkpoint_files/"
            "gnomad.v4.1.1.tmp.ht"
        )


class TestFilterForTest:
    """Test the filter_for_test function."""

    @staticmethod
    def _ht(loci: List[str], build: str) -> hl.Table:
        """Return a Table keyed by locus with one row per locus string."""
        return hl.Table.parallelize(
            [{"locus": hl.parse_locus(l, reference_genome=build)} for l in loci],
            hl.tstruct(locus=hl.tlocus(build)),
            key="locus",
        )

    @pytest.mark.parametrize(
        "build, keep, drop",
        [
            (
                "GRCh38",
                ["chr20:100", "chrX:100", "chrY:2800000"],
                ["chr1:100", "chr19:100", "chr21:100"],
            ),
            ("GRCh37", ["20:100", "X:100", "Y:2800000"], ["1:100", "21:100"]),
        ],
    )
    def test_default_keeps_chr20_x_y(
        self, build: str, keep: List[str], drop: List[str]
    ) -> None:
        """Test that the default test subset is chr20, chrX, and chrY."""
        ht = filter_for_test(self._ht(keep + drop, build))
        assert sorted(str(r.locus) for r in ht.collect()) == sorted(keep)

    def test_gene_list(self) -> None:
        """Test that the gene list keeps only loci inside the test genes.

        The gene intervals are half-open, so the listed end position is excluded.
        """
        keep = ["chr1:55039447", "chr20:50931436", "chrX:13740000", "chrY:2950000"]
        drop = ["chr1:55039446", "chr20:50931437", "chr20:100", "chrX:100"]
        ht = filter_for_test(self._ht(keep + drop, "GRCh38"), use_gene_list=True)
        assert sorted(str(r.locus) for r in ht.collect()) == sorted(keep)

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "filter_for_test compares the ReferenceGenome object from"
            " get_reference_genome to the string 'GRCh37', which is always False,"
            " so GRCh37 Tables get the GRCh38 gene intervals and"
            " parse_locus_interval fails on 'chr1:...'."
        ),
    )
    def test_gene_list_grch37(self) -> None:
        """Test that the gene list uses the GRCh37 gene intervals for GRCh37 Tables."""
        keep = ["1:55505149", "X:13760000"]
        drop = ["1:55505148", "20:100"]
        ht = filter_for_test(self._ht(keep + drop, "GRCh37"), use_gene_list=True)
        assert sorted(str(r.locus) for r in ht.collect()) == sorted(keep)

    def test_keeps_non_locus_fields(self) -> None:
        """Test that filtering keeps the Table's other fields and key."""
        ht = self._ht(["chr20:100", "chr1:100"], "GRCh38")
        ht = filter_for_test(ht.annotate(value=hl.int32(1)))
        assert list(ht.key) == ["locus"]
        assert [r.value for r in ht.collect()] == [1]
