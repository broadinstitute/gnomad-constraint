import hail as hl
from gnomad.resources.grch38.reference_data import clinvar
from gnomad.utils.constraint import oe_confidence_interval
from hail.utils.misc import divide_null

from gnomad_constraint.experimental.promis3d.resources import (
    CURRENT_VERSION,
    get_all_rmc_ht,
    get_constraint_metrics_ht,
    get_context_preprocessed_ht,
    get_cosmis_score_ht,
    get_cosmis_score_tsv,
    get_interpro_annotations,
    get_kaplanis_sig_variants,
    get_kaplanis_variants,
    get_kaplanis_variants_annotated_ht,
)


def process_cosmis_tsv_for_model(
    model: str, version: str = CURRENT_VERSION
) -> hl.Table:
    """
    Import, process, and write the COSMIS TSV file for a specified structure model as a Hail Table.

    This function:
    1. Uses the `get_cosmis_score_tsv` function to locate the COSMIS `.tsv.gz` file.
    2. Imports the TSV with `force=True` and `impute=True`.
    3. Renames/transforms key columns:
       - `transcript_id` from `enst_id`
       - `residue_index` from `uniprot_pos - 1` (to use 0-based indexing)
       - Ensures `cossyn` is a float
    4. Keys the table by (`uniprot_id`, `transcript_id`, `residue_index`).
    5. Writes the resulting Hail Table using the `get_cosmis_score_ht` resource path with checkpointing.

    :param model: Structure model source to process ('alphafold', 'swiss_model', or 'pdb').
    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: The processed and checkpointed Hail Table.
    """
    input_path = get_cosmis_score_tsv(model, version)
    output_path = get_cosmis_score_ht(model, version).path

    ht = hl.import_table(input_path, force=True, impute=True)
    ht = ht.transmute(
        transcript_id=ht.enst_id,
        residue_index=ht.uniprot_pos - 1,
        cossyn=hl.float(ht.cossyn),
    )
    ht = ht.key_by("uniprot_id", "transcript_id", "residue_index")
    ht = ht.checkpoint(output_path, _read_if_exists=True)

    return ht


def process_all_rmc_ht(version: str = CURRENT_VERSION) -> hl.Table:
    """
    Read, annotate, and write the RMC Hail Table with updated p-values and per-locus expansion.

    This function:
    1. Reads the RMC Hail Table using `get_all_rmc_ht`.
    2. Annotates:
       - `section_oe_upper`: Upper bound of the observed/expected confidence interval.
       - `section_p_value`: Chi-square p-value.
       - `locus`: List of loci covered by the interval.
    3. Explodes the `locus` array into multiple rows.
    4. Keys the table by `locus` and `transcript`.
    5. Writes back to the same location using `checkpoint`.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: The processed and checkpointed Hail Table.
    """
    rmc_ht = get_all_rmc_ht(version).ht()

    rmc_ht = rmc_ht.annotate(
        section_oe=divide_null(rmc_ht.section_obs, rmc_ht.section_exp),
        section_oe_upper=(
            hl.qchisqtail(1 - 0.05 / 2, 2 * (rmc_ht.section_obs + 1), lower_tail=True)
            / (2 * rmc_ht.section_exp)
        ),
        section_oe_ci=oe_confidence_interval(rmc_ht.section_obs, rmc_ht.section_exp),
        section_p_value=hl.pchisqtail(rmc_ht.section_chisq, 1),
        locus=hl.range(
            rmc_ht.interval.start.position, rmc_ht.interval.end.position + 1
        ).map(
            lambda x: hl.locus(
                rmc_ht.interval.start.contig, x, reference_genome="GRCh38"
            )
        ),
    )

    rmc_ht = rmc_ht.explode("locus").key_by("locus", "transcript")

    return rmc_ht.checkpoint(
        # hl.utils.new_temp_file("all_rmc_ht")
        "gs://gnomad-tmp-4day/all_rmc_ht-iBnEVygR5ID7isSXAZkUgW",
        _read_if_exists=True,
        # overwrite=True
    )


def process_context_ht(ht, version: str = CURRENT_VERSION) -> hl.Table:
    """
    Read, re-key, and restructure the context preprocessed Hail Table.

    This function:
    1. Reads the context table from `get_context_preprocessed_ht`.
    2. Extracts population labels from the `exomes_freq_meta` field.
    3. Re-keys the table by (`locus`, `alleles`).
    4. Selects and renames key fields:
       - Coverage statistics (`mean`, `median_approx`, `AN`, `percent_AN`)
       - Calibrated mutation frequencies per population

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: The processed and checkpointed Hail Table.
    """
    context_ht = get_context_preprocessed_ht(version).ht()
    ht = context_ht.semi_join(ht.select().key_by("locus", "alleles"))

    pops = hl.eval(ht.exomes_freq_meta.map(lambda x: x.get("gen_anc", "total")))

    ht = ht.select(
        "context",
        "ref",
        "alt",
        "was_flipped",
        "transition",
        "cpg",
        "mutation_type",
        "methylation_level",
        "vep",
        exomes_coverage=hl.struct(
            mean=ht.coverage.exomes.mean,
            median_approx=ht.coverage.exomes.median_approx,
            AN=ht.AN.exomes,
            percent_AN=ht.exomes_coverage,
        ),
        exomes_freq=hl.struct(
            **{pop: ht.calibrate_mu.exomes_freq[i] for i, pop in enumerate(pops)}
        ),
    )

    return ht.checkpoint(
        # hl.utils.new_temp_file("context_preprocessed_ht"),
        "gs://gnomad-tmp-4day/context_preprocessed_ht-c6rnPg63qfL01rmNYBYnOo",
        _read_if_exists=True,
        # overwrite=True
    )


def process_constraint_metrics_ht(version: str = CURRENT_VERSION) -> hl.Table:
    """
    Read, filter, and restructure the transcript constraint metrics Hail Table.

    This function:
    1. Reads the constraint metrics HT from `get_constraint_metrics_ht`.
    2. Filters to rows with transcript IDs starting with "ENST".
    3. Selects three summary structs (`syn`, `mis`, `lof`) from `constraint_groups`,
       each augmented with its respective `oe_info`.
    4. Keys the table by `transcript`.
    5. Writes the result back using `checkpoint`.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: The processed and checkpointed Hail Table.
    """
    ht = get_constraint_metrics_ht(version).ht()
    ht = ht.filter(ht.transcript.startswith("ENST"))

    ht = ht.select(
        syn=ht.constraint_groups[0]
        .annotate(**ht.constraint_groups[0].oe_info[0])
        .drop("oe_info"),
        mis=ht.constraint_groups[1]
        .annotate(**ht.constraint_groups[1].oe_info[1])
        .drop("oe_info"),
        lof=ht.constraint_groups[5]
        .annotate(**ht.constraint_groups[5].oe_info[5])
        .drop("oe_info"),
    )

    ht = ht.key_by("transcript")

    return ht


def process_kaplanis_variants_ht() -> hl.Table:
    """
    Process Kaplanis de novo missense variants file and lift over loci to GRCh38.

    This function:
    1. Imports the Kaplanis annotated variants TSV using `get_kaplanis_variants`.
    2. Filters to `missense_variant` consequences.
    3. Parses locus and alleles from the `Variant` field.
    4. Selects relevant fields and performs liftover from GRCh37 to GRCh38.
    5. Filters variants to those within significant consensus genes (from `get_kaplanis_sig_variants`).
    6. Keys by (`locus`, `alleles`, `transcript_id`) on GRCh38.
    7. Saves the result to the same path using `checkpoint`.

    :return: The processed and checkpointed Hail Table.
    """
    ht = hl.import_table(get_kaplanis_variants(), min_partitions=10, impute=True)
    ht = ht.filter(ht.Consequence_Kaplanis == "missense_variant")

    variant_expr = ht.Variant.split(":")
    ht = ht.select(
        locus=hl.locus(
            variant_expr[0], hl.int(variant_expr[1]), reference_genome="GRCh37"
        ),
        transcript_id=ht.Transcript_ID,
        gene_id=ht.Gene_ID,
        gene_name=ht.Gene,
        alleles=[variant_expr[2], variant_expr[3]],
        mpc=ht.MPC,
        am_pathogenicity=ht.am_pathogenicity,
        mpc_v2=ht.MPC_v2,
        mpc_v2_Outlier=ht.MPC_v2_Outlier,
        polyphen=ht.PolyPhen,
        sift=ht.Sift,
    )

    # Load significant gene set
    sig_gene_ht = hl.import_table(get_kaplanis_sig_variants())
    sig_genes = sig_gene_ht.aggregate(
        hl.agg.filter(
            (sig_gene_ht.significant == "TRUE")
            | (sig_gene_ht.diagnostic_category == "consensus"),
            hl.agg.collect_as_set(sig_gene_ht.symbol),
        )
    )
    ht = ht.filter(hl.literal(sig_genes).contains(ht.gene_name))

    # Perform liftover
    rg37 = hl.get_reference("GRCh37")
    rg38 = hl.get_reference("GRCh38")
    rg37.add_liftover(
        "gs://hail-common/references/grch37_to_grch38.over.chain.gz", rg38
    )

    ht = ht.annotate(
        grch37_locus=ht.locus,
        lift_locus=hl.liftover(ht.locus, "GRCh38"),
    )

    ht = ht.key_by(
        locus=ht.lift_locus,
        alleles=ht.alleles,
        # transcript_id=ht.transcript_id,
    )

    ht = ht.select(
        "grch37_locus",
        "mpc",
        "am_pathogenicity",
        "mpc_v2",
        "mpc_v2_Outlier",
        "polyphen",
        "sift",
        "gene_name",
    )

    ht = ht.checkpoint(
        get_kaplanis_variants_annotated_ht().path,
        # _read_if_exists=True,
        overwrite=True,
    )

    return ht


def process_interpro_ht(version: str = CURRENT_VERSION) -> hl.Table:
    """
    Process the InterPro annotation file into a Hail Table with per-residue indexing.

    This function:
    1. Imports the InterPro TSV using `get_interpro_annotations`.
    2. Parses start and end residue indices, handles missing values.
    3. Explodes into individual residue indices (0-based).
    4. Selects and renames relevant fields.
    5. Keys by (`uniprot_id`, `transcript_id`, `residue_index`).
    6. Writes the resulting table to the path returned by `get_interpro_ht`.

    :param version: gnomAD version to use. Default is `CURRENT_VERSION`.
    :return: The processed and checkpointed Hail Table.
    """
    ht = hl.import_table(
        get_interpro_annotations(version), impute=True, min_partitions=500
    )

    ht = ht.select(
        uniprot_id=ht["UniProtKB/Swiss-Prot ID"],
        transcript_id=ht["Transcript stable ID"],
        residue_index=hl.range(
            hl.int(hl.or_missing(ht["Interpro start"] != "", ht["Interpro start"])) - 1,
            hl.int(hl.or_missing(ht["Interpro end"] != "", ht["Interpro end"])),
        ),
        interpro_id=ht["Interpro ID"],
        interpro_short_description=ht["Interpro Short Description"],
        interpro_description=ht["Interpro Description"],
    )

    ht = ht.explode("residue_index")
    ht = ht.key_by("uniprot_id", "transcript_id", "residue_index")

    return ht


def process_clinvar_ht() -> hl.Table:
    """
    Process ClinVar HT by extracting info fields, filtering to missense variants, and keying by gene.

    This function:
    1. Reads ClinVar GRCh38 Hail Table from the public Google Cloud path.
    2. Selects all fields from `info` and the top-level `rsid`.
    3. Filters to variants with at least one consequence marked as `missense_variant` in the MC field.
    4. Keys the table by (`locus`, `alleles`, `gene`), where `gene` is extracted from the `GENEINFO` field.
    5. Saves the result using checkpoint to the path from `get_clinvar_missense_ht`.

    :return: Processed and checkpointed Hail Table.
    """
    ht = clinvar.ht()
    ht = ht.select("rsid", **ht.info)
    ht = ht.filter(hl.any(ht.MC.map(lambda x: x.split("\\|")[1] == "missense_variant")))
    ht = ht.key_by("locus", "alleles", gene=ht.GENEINFO.split(":")[0])

    return ht.checkpoint(
        # hl.utils.new_temp_file("clinvar_missense_ht")
        "gs://gnomad-tmp-4day/clinvar_missense_ht-mKswH0ptCtM4PsuoQOUkCw",
        _read_if_exists=True,
        # overwrite=True
    )
