"""Tests for the resource utility module."""

from typing import Any, Callable, Dict, List, Union

import gnomad.resources.grch38.gnomad as gnomad_grch38
import gnomad.resources.grch38.reference_data as ref_grch38
import hail as hl
import pytest
from gnomad.resources.resource_utils import ExpressionResource, TableResource
from gnomad_qc.resource_utils import PipelineResourceCollection

from gnomad_constraint.resources.constants import (
    GENCODE_VERSION_MAP,
    VEP_VERSION_MAP,
    VERSIONS,
)
from gnomad_constraint.resources.resource_utils import (
    check_param_scope,
    filter_for_test,
    get_adj_r_ht,
    get_aggregated_expected,
    get_aggregated_per_variant_expected,
    get_annotated_context_ht,
    get_constraint_data,
    get_constraint_group_ht,
    get_constraint_metrics_dataset,
    get_constraint_resources,
    get_constraint_root,
    get_coverage_ht,
    get_dnm_ht,
    get_gene_quality_metrics_ht,
    get_genomic_features_ht,
    get_lof_threshold_tsv_path,
    get_logging_path,
    get_models,
    get_mutation_ht,
    get_per_variant_expected_dataset,
    get_pre_rank_constraint_metrics_dataset,
    get_preprocessed_ht,
    get_regional_rr_ht,
    get_release_constraint_ht,
    get_release_constraint_tsv_path,
    get_release_downsampling_tsv_path,
    get_release_mutation_ht,
    get_release_mutation_tsv_path,
    get_rr_gene_flags_ht,
    get_rr_models_ht,
    get_rr_penalty_tuning_ht,
    get_rr_site_features_ht,
    get_rr_training_ht,
    get_sites_resource,
    get_syn_adj_r_ht,
    get_training_dataset,
    get_training_tsv_path,
    get_vep_context_ht,
)

PROD_ROOT = "gs://gnomad/v4.1.1/constraint"
TEST_ROOT = "gs://gnomad-tmp/gnomad_v4.1.1_testing/constraint"


def _path(resource: Union[str, TableResource]) -> str:
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
            ({"extension": "csv"}, "extension must be one of"),
            ({"version": "3.1"}, "The requested version doesn't exist"),
            ({"version": "2.1.1"}, "The requested version doesn't exist"),
        ],
    )
    def test_invalid_argument(self, kwargs: Dict[str, str], message: str) -> None:
        """Test that each out-of-scope argument raises an error naming it."""
        with pytest.raises(ValueError, match=message):
            check_param_scope(**kwargs)

    def test_valid_arguments(self) -> None:
        """Test that valid arguments return None."""
        assert (
            check_param_scope(
                version="4.1.1",
                model_type="plateau",
                extension="ht",
                data_type="exomes",
            )
            is None
        )


class TestGetVepContextHt:
    """Test the get_vep_context_ht function."""

    @pytest.mark.parametrize("version", ["4.0", "4.1.1"])
    def test_version_mapping(self, version: str) -> None:
        """Test that v4 versions use the VEP 105 context."""
        expected = ref_grch38.vep_context.versions["105"]
        assert get_vep_context_ht(version).path == expected.path

    @pytest.mark.parametrize("version", ["2.1.1", "3.1"])
    def test_invalid_version(self, version: str) -> None:
        """Test that a version outside v4 raises an error."""
        with pytest.raises(ValueError, match="The requested version doesn't exist"):
            get_vep_context_ht(version)


class TestGencodeVersion:
    """Test the GENCODE version used for each pipeline version."""

    @pytest.mark.parametrize("version", VERSIONS)
    def test_vep_version_has_gencode_resource(self, version: str) -> None:
        """
        Test that each version's VEP version maps to a GENCODE resource.

        ``get_gencode_ht`` takes the GENCODE version from the context Table's VEP
        version, not from gnomad_methods' default GENCODE version, which can change
        independently of the VEP annotations.
        """
        vep_version = VEP_VERSION_MAP[version]
        assert (
            get_vep_context_ht(version).path
            == ref_grch38.vep_context.versions[vep_version].path
        )
        assert GENCODE_VERSION_MAP[vep_version] in ref_grch38.gencode.versions


class TestGetSitesResource:
    """Test the get_sites_resource function."""

    @pytest.mark.parametrize(
        "data_type, version, expected",
        [
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
            ("exomes", "4.1.1", gnomad_grch38.coverage("exomes")),
            ("genomes", "4.1.1", gnomad_grch38.coverage("genomes")),
        ],
    )
    def test_build_mapping(self, data_type: str, version: str, expected: Any) -> None:
        """Test that the coverage resource is gnomAD's GRCh38 coverage."""
        resource = get_coverage_ht(data_type, version)
        assert resource.path == expected.path
        assert list(resource.versions) == list(expected.versions)


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
            (get_dnm_ht, "imported_data/gnomad.v4.1.1.dnm.ht"),
            (
                get_genomic_features_ht,
                "imported_data/gnomad.v4.1.1.genomic_features13_1kb.ht",
            ),
            (get_regional_rr_ht, "regional_rate/gnomad.v4.1.1.syn_dnm_rr_1kb.ht"),
            (
                get_rr_site_features_ht,
                "regional_rate/gnomad.v4.1.1.rr_site_features.ht",
            ),
            (get_rr_training_ht, "regional_rate/gnomad.v4.1.1.rr_training.ht"),
            (
                get_rr_penalty_tuning_ht,
                "regional_rate/gnomad.v4.1.1.rr_penalty_tuning.ht",
            ),
            (get_rr_models_ht, "regional_rate/gnomad.v4.1.1.rr_models.ht"),
        ],
    )
    def test_path(self, func: Callable, path: str) -> None:
        """Test the production and test paths of each output."""
        assert _path(func(version="4.1.1")) == f"{PROD_ROOT}/{path}"
        assert _path(func(version="4.1.1", test=True)) == f"{TEST_ROOT}/{path}"

    @pytest.mark.parametrize("func, sub_dir, name", VEP_FUNCS)
    def test_transcript_consequences_subdir(
        self, func: Callable, sub_dir: str, name: str
    ) -> None:
        """Test that the outputs are under ``transcript_consequences`` and post-fixes apply.

        The subdirectory is left from when the VEP annotation was an option, and is
        kept so the paths of existing outputs don't change.
        """
        out = func(version="4.1.1", directory_post_fix="ukb", path_post_fix="v2")
        assert out.path == (
            f"gs://gnomad/v4.1.1/constraint_ukb/{sub_dir}/transcript_consequences/"
            f"gnomad.v4.1.1.{name}.v2.ht"
        )

    @pytest.mark.parametrize("func, sub_dir, name", VEP_FUNCS)
    def test_vep_annotation_is_not_an_option(
        self, func: Callable, sub_dir: str, name: str
    ) -> None:
        """Test that a VEP annotation can no longer be chosen.

        ``worst_csq_by_gene`` never ran end to end on v4: the metrics, ranking, and
        release steps need transcripts. See issue #62.
        """
        with pytest.raises(TypeError):
            func(custom_vep_annotation="worst_csq_by_gene", version="4.1.1")


class TestGetLoggingPath:
    """Test the get_logging_path function."""

    def test_always_under_test_root(self) -> None:
        """Test that log files always go under the test root."""
        assert get_logging_path("constraint", version="4.1.1") == (
            f"{TEST_ROOT}/logging/gnomad.v4.1.1.constraint.log"
        )

    def test_post_fixes(self) -> None:
        """Test that the post-fixes keep the logs of separate runs apart."""
        path = get_logging_path(
            "constraint", version="4.1.1", directory_post_fix="x", path_post_fix="v2"
        )
        assert path == (
            "gs://gnomad-tmp/gnomad_v4.1.1_testing/constraint_x/logging/"
            "gnomad.v4.1.1.constraint.v2.log"
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

    def test_default_keeps_chr20_x_y(self) -> None:
        """Test that the default test subset is chr20, chrX, and chrY."""
        keep = ["chr20:100", "chrX:100", "chrY:2800000"]
        drop = ["chr1:100", "chr19:100", "chr21:100"]
        ht = filter_for_test(self._ht(keep + drop, "GRCh38"))
        assert sorted(str(r.locus) for r in ht.collect()) == sorted(keep)

    def test_gene_list(self) -> None:
        """Test that the gene list keeps only loci inside the test genes.

        The gene intervals are half-open, so the listed end position is excluded.
        """
        keep = ["chr1:55039447", "chr20:50931436", "chrX:13740000", "chrY:2950000"]
        drop = ["chr1:55039446", "chr20:50931437", "chr20:100", "chrX:100"]
        ht = filter_for_test(self._ht(keep + drop, "GRCh38"), use_gene_list=True)
        assert sorted(str(r.locus) for r in ht.collect()) == sorted(keep)

    def test_keeps_non_locus_fields(self) -> None:
        """Test that filtering keeps the Table's other fields and key."""
        ht = self._ht(["chr20:100", "chr1:100"], "GRCh38")
        ht = filter_for_test(ht.annotate(value=hl.int32(1)))
        assert list(ht.key) == ["locus"]
        assert [r.value for r in ht.collect()] == [1]


class TestGetConstraintResources:
    """Test the pipeline steps built by the get_constraint_resources function."""

    @staticmethod
    def _get_resources(
        skip_pre_rank_metrics: bool = False,
    ) -> PipelineResourceCollection:
        """Return the constraint pipeline resources for a test run of v4.1.1."""
        return get_constraint_resources(
            version="4.1.1",
            overwrite=False,
            test=True,
            skip_pre_rank_metrics=skip_pre_rank_metrics,
        )

    @pytest.fixture
    def resources(self) -> PipelineResourceCollection:
        """Return the constraint pipeline resources for a test run of v4.1.1."""
        return self._get_resources()

    def test_export_steps_check_only_their_own_outputs(
        self, resources: PipelineResourceCollection
    ) -> None:
        """
        Test that each release export flag checks only the files it writes.

        Regression test: both export flags used to share one step, so once either
        had run, the other refused to run without ``--overwrite`` because the first
        one's files already existed.
        """
        assert resources.export_release_tsv.output_resources == {
            "--export-release-tsv": [
                get_release_constraint_tsv_path(version="4.1.1", test=True),
                get_lof_threshold_tsv_path(version="4.1.1", test=True),
            ]
        }
        assert resources.export_release_downsampling_tsv.output_resources == {
            "--export-release-downsampling-tsv": [
                get_release_downsampling_tsv_path(version="4.1.1", test=True)
            ]
        }

    @pytest.mark.parametrize("skip_pre_rank_metrics", [False, True])
    def test_inputs_are_keyed_by_the_flag_that_creates_them(
        self, skip_pre_rank_metrics: bool
    ) -> None:
        """
        Test that every step's inputs are keyed by the flag that creates them.

        When an input is missing, gnomad_qc's error tells the user to add the key to
        the command line, so the key has to be a pipeline flag. Regression test: the
        preprocessing step was named "preprocess data for downstream steps", and
        with ``--skip-pre-rank-metrics`` the pre-rank Table was keyed by
        "pre-rank constraint metrics". The gnomAD inputs, the published files the
        rr imports read, and the constraint metrics Table the rr site features take
        gene flags from aren't created by the pipeline, so they are exempt.
        """
        resources = self._get_resources(skip_pre_rank_metrics)
        keys = {
            key
            for step in resources.pipeline_steps.values()
            for key in step.input_resources
        }
        assert "--preprocess-data" in keys
        external = {
            "gnomAD resources",
            "published DNMs",
            "published genomic features",
            "constraint metrics with gene flags",
        }
        assert all(key.startswith("--") for key in keys - external)
        if skip_pre_rank_metrics:
            assert list(resources.compute_constraint_metrics.input_resources) == [
                "--compute-constraint-metrics"
            ]

    def test_test_run_writes_only_under_the_test_root(
        self, resources: PipelineResourceCollection
    ) -> None:
        """
        Test that every step's outputs in a test run are under the test root.

        Regression test for 3a5c64b: the release and gene quality outputs ignored
        ``test``, so a ``--test --overwrite`` run wrote over the production files.
        Checking the resources built by ``get_constraint_resources`` also catches
        it no longer passing ``test`` to one of the path functions.
        """
        paths = [
            _path(resource)
            for step in resources.pipeline_steps.values()
            for outputs in step.output_resources.values()
            for resource in outputs
        ]
        assert paths
        assert [p for p in paths if not p.startswith(f"{TEST_ROOT}/")] == []

    @pytest.mark.parametrize("use_aggregated_expected", [False, True])
    def test_constraint_group_input(self, use_aggregated_expected: bool) -> None:
        """
        Test that --aggregate-by-constraint-groups declares the input it reads.

        Regression test: with ``--use-aggregated-expected``, the step read the
        ``--apply-models-aggregated`` output without declaring it, so its resources
        weren't checked.
        """
        resources = get_constraint_resources(
            version="4.1.1",
            overwrite=False,
            test=True,
            use_aggregated_expected=use_aggregated_expected,
        )
        if use_aggregated_expected:
            step = "--apply-models-aggregated"
            path = get_aggregated_expected(version="4.1.1", test=True).path
        else:
            step = "--aggregate-per-variant-expected"
            path = get_aggregated_per_variant_expected(version="4.1.1", test=True).path
        inputs = resources.aggregate_by_constraint_groups.input_resources
        assert list(inputs) == [step]
        assert [_path(r) for r in inputs[step]] == [path]

    def test_inputs_read_by_steps_are_declared(
        self, resources: PipelineResourceCollection
    ) -> None:
        """
        Test that the Tables steps read are declared inputs, so they're checked.

        Regression test: ``--prepare-context-ht`` read the adj_r and syn_adj_r Tables,
        and ``--compute-gene-quality-metrics`` the segdup and LCR intervals, without
        declaring them, so a missing one failed on read instead of in the resource
        check.
        """

        def _input_paths(step: Any) -> set:
            return {_path(r) for rs in step.input_resources.values() for r in rs}

        assert {get_adj_r_ht().path, get_syn_adj_r_ht().path} <= _input_paths(
            resources.prepare_context
        )
        assert {
            ref_grch38.seg_dup_intervals.path,
            ref_grch38.lcr_intervals.path,
        } <= _input_paths(resources.compute_gene_quality_metrics)


class TestGetConstraintResourcesImportSteps:
    """Test the rr import steps built by the get_constraint_resources function."""

    @pytest.mark.parametrize("test", [False, True])
    def test_outputs(self, test: bool) -> None:
        """
        Test that each step writes its Table, under the test root in a test run.

        The import outputs take neither post-fix.
        """
        resources = get_constraint_resources(
            "4.1.1",
            overwrite=False,
            test=test,
            directory_post_fix="d",
            path_post_fix="p",
        )
        assert (
            resources.import_dnms.dnm_ht.path
            == get_dnm_ht(version="4.1.1", test=test).path
        )
        assert resources.import_genomic_features.genomic_features_ht.path == (
            get_genomic_features_ht(version="4.1.1", test=test).path
        )

    def test_inputs(self) -> None:
        """Test that the inputs the import steps read are the published files."""
        resources = get_constraint_resources("4.1.1", overwrite=False, test=True)
        dnm_dir = "gs://regional_missense_constraint/resources/GRCh38/reference_data"
        assert resources.import_dnms.asd_dnm_path == (
            f"{dnm_dir}/Satterstrom_2026_ASD_dnv_SNV_indel.txt.gz"
        )
        assert resources.import_dnms.ddd_dnm_path == (
            f"{dnm_dir}/Satterstrom_2026_DDD_dnv_SNV_indel.txt.gz"
        )
        assert resources.import_genomic_features.published_features_path == (
            "gs://gnomad-nc-constraint-v31-paper/misc/genomic_features13_genome_1kb.txt"
        )


class TestGetConstraintResourcesRrSteps:
    """Test the rr steps built by the get_constraint_resources function."""

    def test_inputs(self) -> None:
        """
        Test that the rr steps read the imports, each other's outputs, and the flags.

        The gene flags come from the v4.1.2 constraint metrics Table whatever the
        version. The training set is made from the site features, and the tuning and
        the fit read it.
        """
        resources = get_constraint_resources("4.1.1", overwrite=False, test=True)
        sites = resources.create_rr_site_features
        assert sites.dnm_ht.path == resources.import_dnms.dnm_ht.path
        for step in [sites, resources.compute_rr]:
            assert (
                step.genomic_features_ht.path
                == resources.import_genomic_features.genomic_features_ht.path
            )
        assert sites.gene_flags_ht.path == get_rr_gene_flags_ht().path
        assert "/v4.1.2/" in get_rr_gene_flags_ht().path
        training = resources.create_rr_training_set
        assert training.rr_site_features_ht.path == sites.rr_site_features_ht.path
        for step in [resources.tune_rr_penalty, resources.train_rr_models]:
            assert step.rr_training_ht.path == training.rr_training_ht.path
        train = resources.train_rr_models
        assert (
            train.rr_penalty_tuning_ht.path
            == resources.tune_rr_penalty.rr_penalty_tuning_ht.path
        )
        assert (
            resources.compute_rr.rr_models_ht.path
            == resources.train_rr_models.rr_models_ht.path
        )

    @pytest.mark.parametrize("test", [False, True])
    @pytest.mark.parametrize("directory_post_fix", [None, "d"])
    @pytest.mark.parametrize("path_post_fix", [None, "p"])
    def test_inputs_match_the_constraint_steps(
        self, test: bool, directory_post_fix: str, path_post_fix: str
    ) -> None:
        """
        Test that the rr steps read the Tables of the same constraint pipeline run.

        For the same version, test flag, and post-fixes, the site features come from
        the annotated context Table that ``--prepare-context-ht`` writes, and rr is
        normalized on the per-variant expected Table that
        ``--aggregate-per-variant-expected`` aggregates. ``--compute-rr`` writes the
        rr Table to its getter's path with the same settings, which
        ``--apply-regional-rr`` reads.
        """
        kwargs = {
            "test": test,
            "directory_post_fix": directory_post_fix,
            "path_post_fix": path_post_fix,
        }
        resources = get_constraint_resources("4.1.1", overwrite=False, **kwargs)
        assert (
            resources.create_rr_site_features.annotated_context_ht.path
            == resources.prepare_context.annotated_context_ht.path
        )
        assert (
            resources.compute_rr.per_variant_apply_ht.path
            == resources.aggregate_per_variant_expected.per_variant_apply_ht.path
        )
        assert (
            resources.compute_rr.regional_rr_ht.path
            == get_regional_rr_ht(version="4.1.1", **kwargs).path
        )
