"""
Script to compute Site Frequency Spectrum (SFS) tables from constraint pipeline data.

This script generates SFS tables for each combination of:
- Functional categories: synonymous (syn), missense (mis), loss-of-function (lof)
- Mutation type categories: transversion, transition, CpG (methylated and non-methylated)
- Optional: Genetic ancestry groups (afr, amr, eas, nfe, sas)

The output format is tables with:
- Rows: Different downsampling levels (AC_ds10, AC_ds100, etc.)
- Columns: Allele count bins (0, 1, 2, ..., 10, 10-100, 100-1000, 1000-10000, >10000)
- Values: Number of variants with that allele count in that downsampling

When grouping by ancestry, separate tables are generated for each ancestry group.

The output is modeled after the SFS tables sent to Julia Goodrich by Jeremy Guez.
"""

import argparse
import logging
from typing import Dict, List, Optional, Union

import hail as hl
from gnomad.utils.vep import CSQ_ORDER, LOF_CSQ_SET, get_most_severe_consequence_expr

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("compute_sfs_tables")
logger.setLevel(logging.INFO)


def create_ac_group_expr(ac_expr: hl.expr.Int32Expression) -> hl.expr.StringExpression:
    """
    Create expression to assign variants to allele count bins.

    :param ac_expr: Expression for allele count
    :return: Expression that assigns variants to AC bins
    """
    return (
        hl.case()
        .when(ac_expr == 0, "0")
        .when(ac_expr == 1, "1")
        .when(ac_expr == 2, "2")
        .when(ac_expr == 3, "3")
        .when(ac_expr == 4, "4")
        .when(ac_expr == 5, "5")
        .when(ac_expr == 6, "6")
        .when(ac_expr == 7, "7")
        .when(ac_expr == 8, "8")
        .when(ac_expr == 9, "9")
        .when(ac_expr == 10, "10")
        .when((ac_expr > 10) & (ac_expr <= 100), "10-100")
        .when((ac_expr > 100) & (ac_expr <= 1000), "100-1000")
        .when((ac_expr > 1000) & (ac_expr <= 10000), "1000-10000")
        .when(ac_expr > 10000, ">10000")
        .or_missing()
    )


def get_most_severe_consequence_expr(
    csq_expr: hl.expr.ArrayExpression,
    csq_order: Optional[List[str]] = None,
    csq_field: str = "most_severe_consequence",
) -> Union[hl.expr.StringExpression, hl.expr.StructExpression]:
    """
    Get the most severe consequence from a collection of consequences.

    This is for a given transcript, as there are often multiple annotations for a single
    transcript: e.g. splice_region_variant&intron_variant -> splice_region_variant

    :param csq_expr: ArrayExpression of consequences.
    :param csq_order: Optional list indicating the order of VEP consequences, sorted
        from high to low impact. Default is None, which uses the value of the
        `CSQ_ORDER` global.
    :param csq_field: Field name for consequence in the struct. Default is
      "most_severe_consequence".
    :return: Most severe consequence in `csq_expr`.
    """
    if csq_order is None:
        csq_order = CSQ_ORDER
    csqs = hl.literal(csq_order)

    if csq_expr.dtype == hl.tarray(hl.tstr):
        return csqs.find(lambda c: csq_expr.contains(c))
    else:
        ms_csq = csqs.find(lambda c: csq_expr.map(lambda x: x[csq_field]).contains(c))
        return csq_expr.filter(lambda x: x[csq_field] == ms_csq).first()


def preprocess_constraint_data(
    ht: hl.Table, gen_ancs: Optional[List[str]] = None
) -> hl.Table:
    """
    Preprocess constraint data for SFS computation.

    This function:

        - Loads the constraint data
        - Selects relevant annotations
        - Creates AC group expressions
        - Explodes frequency data by downsampling level and ancestry group
        - Saves the preprocessed data

    Compute the allele count of variants in the same way as we do for observed and
    possible variants in the constraint pipeline:

        - observed_variants: This annotation is an array, where each element
          corresponds to whether the variant is observed in the genomes dataset for
          the frequency group at the corresponding index in the `genomes_freq` array.
          Must PASS genome filters, have AC <= 'ac_cutoff' at the specified
          `downsampling_level`, and have a genome mean coverage >= `min_cov` and <=
          `max_cov`. The boolean value is stored as an integer (0 or 1).
        - possible_variants: Whether the variant is considered a possible variant in
          the genomes dataset. This includes variants not in the genome dataset
          (genome AF undefined), or also considered in the observed variant set. The
          boolean value is stored as an integer (0 or 1).

    :param ht: Input constraint table.
    :param gen_ancs: Optional list of genetic ancestry groups to include.
        If None, only "global" is included. Default is None.
    :return: Preprocessed Hail table
    """
    # Define consequences to keep.
    csq_to_keep = ["missense_variant", "synonymous_variant"] + list(LOF_CSQ_SET)

    # Filter to variants with transcript consequences.
    ht = ht.filter(
        hl.is_defined(ht.calibrate_mu.possible_variants)
        & (hl.set(csq_to_keep).contains(ht.vep.most_severe_consequence))
    )

    # Get frequency metadata.
    freq_meta = hl.eval(ht.exomes_freq_meta)

    # Build list of frequency metadata entries to include.
    freq_meta_include = [(0, "AC_ds730947", "global")]

    # Include global downsamplings.
    freq_meta_include += [
        (i, f"AC_ds{m['downsampling']}", "global")
        for i, m in enumerate(freq_meta)
        if "downsampling" in m.keys()
        and "gen_anc" in m.keys()
        and m["gen_anc"] == "global"
    ]

    # Include ancestry-specific downsamplings if requested.
    if gen_ancs is not None:
        freq_meta_include += [
            (i, f"AC_ds{m['downsampling']}", gen_anc)
            for gen_anc in gen_ancs
            for i, m in enumerate(freq_meta)
            if "downsampling" in m.keys()
            and "gen_anc" in m.keys()
            and m["gen_anc"] == gen_anc
        ]

    # Annotate transcript consequences.
    ht = ht.annotate(
        transcript_consequences=ht.vep.transcript_consequences.map(
            lambda x: x.select(
                "transcript_id",
                "gene_id",
                "gene_symbol",
                "biotype",
                "most_severe_consequence",
                "mane_select",
                "canonical",
                "lof",
            )
        ).filter(
            lambda x: (
                x.transcript_id.startswith("ENST")
                & (x.biotype == "protein_coding")
                & hl.set(csq_to_keep).contains(x.most_severe_consequence)
            )
        ),
        freq=[
            hl.struct(
                ac_group=create_ac_group_expr(
                    hl.if_else(
                        hl.is_defined(ht.filters.exomes)
                        & (ht.filters.exomes.length() == 0),
                        ht.calibrate_mu.exomes_freq[i].AC,
                        0,
                    )
                ),
                downsampling=m,
                gen_anc=gen_anc,
            )
            for i, m, gen_anc in freq_meta_include
        ],
    )

    # Filter to variants with transcript consequences.
    ht = ht.filter(ht.transcript_consequences.length() > 0)

    # Create downsampling annotations.
    logger.info("Creating downsampling annotations")
    csq_expr = ht.transcript_consequences
    ht = ht.select(
        "transition",
        "cpg",
        "mutation_type",
        "methylation_level",
        "freq",
        most_severe_consequence=get_most_severe_consequence_expr(csq_expr),
        canonical_most_severe_consequence=get_most_severe_consequence_expr(
            csq_expr.filter(lambda x: hl.or_else(x.canonical == 1, False))
        ),
        mane_select_most_severe_consequence=get_most_severe_consequence_expr(
            csq_expr.filter(lambda x: hl.is_defined(x.mane_select))
        ),
    )

    # Explode frequency data.
    logger.info("Exploding frequency data by downsampling level")
    ht = ht.explode(ht.freq)
    ht = ht.annotate(**ht.freq).drop("freq")

    return ht


def aggregate_sfs_data(ht: hl.Table, group_by_ancestry: bool = False) -> hl.Table:
    """
    Aggregate SFS data by downsampling level, AC group, and optionally ancestry group.

    :param ht: Preprocessed Hail table.
    :param group_by_ancestry: Whether to group by genetic ancestry. Default is False.
    :return: Aggregated Hail table.
    """
    group_keys = ["downsampling", "ac_group"]
    if group_by_ancestry:
        logger.info(
            "Aggregating SFS data by downsampling level, AC group, and ancestry group"
        )
        group_keys += ["gen_anc"]
    else:
        logger.info("Aggregating SFS data by downsampling level and AC group")

    agg_ht = (
        ht.group_by(*group_keys)
        .aggregate(
            n=hl.struct(
                total=hl.agg.count(),
                most_severe_consequence=hl.agg.group_by(
                    ht.most_severe_consequence.select(
                        "most_severe_consequence",
                        "lof",
                        mutation_type=ht.mutation_type,
                        transition=ht.transition,
                        cpg=ht.cpg,
                        methylation_level=ht.methylation_level,
                    ),
                    hl.agg.count(),
                ),
                canonical_most_severe_consequence=hl.agg.group_by(
                    ht.canonical_most_severe_consequence.select(
                        "most_severe_consequence",
                        "lof",
                        mutation_type=ht.mutation_type,
                        transition=ht.transition,
                        cpg=ht.cpg,
                        methylation_level=ht.methylation_level,
                    ),
                    hl.agg.count(),
                ),
                mane_select_most_severe_consequence=hl.agg.group_by(
                    ht.mane_select_most_severe_consequence.select(
                        "most_severe_consequence",
                        "lof",
                        mutation_type=ht.mutation_type,
                        transition=ht.transition,
                        cpg=ht.cpg,
                        methylation_level=ht.methylation_level,
                    ),
                    hl.agg.count(),
                ),
            )
        )
        .cache()
        .naive_coalesce(10)
    )

    return agg_ht


def get_sfs_grouping_definitions() -> Dict[str, hl.StructExpression]:
    """
    Get definitions for SFS table groupings.

    Returns a dictionary mapping output filenames to their filtering criteria.

    :return: Dictionary of filename -> grouping criteria
    """
    lof_groups = ["HC", "LC", hl.missing(hl.tstr)]

    return {
        "SFS_mis_nonCpGtransition.txt.gz": hl.Struct(
            most_severe_consequence=["missense_variant"],
            lof=[hl.missing(hl.tstr)],
            mutation_type=["non-CpG transition"],
            transition=[True],
            cpg=[False],
            methylation_level=[0, 1, 2],
        ),
        "SFS_mis_transversion.txt.gz": hl.Struct(
            most_severe_consequence=["missense_variant"],
            lof=[hl.missing(hl.tstr)],
            mutation_type=["transversion"],
            transition=[False],
            cpg=[False],
            methylation_level=[0, 1, 2],
        ),
        "SFS_mis_methylCpG.txt.gz": hl.Struct(
            most_severe_consequence=["missense_variant"],
            lof=[hl.missing(hl.tstr)],
            mutation_type=["CpG"],
            transition=[True],
            cpg=[True],
            methylation_level=[2],
        ),
        "SFS_mis_nonmethylCpG.txt.gz": hl.Struct(
            most_severe_consequence=["missense_variant"],
            lof=[hl.missing(hl.tstr)],
            mutation_type=["CpG"],
            transition=[True],
            cpg=[True],
            methylation_level=[0],
        ),
        "SFS_syn_nonCpGtransition.txt.gz": hl.Struct(
            most_severe_consequence=["synonymous_variant"],
            lof=[hl.missing(hl.tstr)],
            mutation_type=["non-CpG transition"],
            transition=[True],
            cpg=[False],
            methylation_level=[0, 1, 2],
        ),
        "SFS_syn_transversion.txt.gz": hl.Struct(
            most_severe_consequence=["synonymous_variant"],
            lof=[hl.missing(hl.tstr)],
            mutation_type=["transversion"],
            transition=[False],
            cpg=[False],
            methylation_level=[0, 1, 2],
        ),
        "SFS_syn_methylCpG.txt.gz": hl.Struct(
            most_severe_consequence=["synonymous_variant"],
            lof=[hl.missing(hl.tstr)],
            mutation_type=["CpG"],
            transition=[True],
            cpg=[True],
            methylation_level=[2],
        ),
        "SFS_syn_nonmethylCpG.txt.gz": hl.Struct(
            most_severe_consequence=["synonymous_variant"],
            lof=[hl.missing(hl.tstr)],
            mutation_type=["CpG"],
            transition=[True],
            cpg=[True],
            methylation_level=[0],
        ),
        "SFS_LoF_nonCpGtransition.txt.gz": hl.Struct(
            most_severe_consequence=LOF_CSQ_SET,
            lof=lof_groups,
            mutation_type=["non-CpG transition"],
            transition=[True],
            cpg=[False],
            methylation_level=[0, 1, 2],
        ),
        "SFS_LoF_transversion.txt.gz": hl.Struct(
            most_severe_consequence=LOF_CSQ_SET,
            lof=lof_groups,
            mutation_type=["transversion"],
            transition=[False],
            cpg=[False],
            methylation_level=[0, 1, 2],
        ),
        "SFS_LoF_methylCpG.txt.gz": hl.Struct(
            most_severe_consequence=LOF_CSQ_SET,
            lof=lof_groups,
            mutation_type=["CpG"],
            transition=[True],
            cpg=[True],
            methylation_level=[2],
        ),
        "SFS_LoF_nonmethylCpG.txt.gz": hl.Struct(
            most_severe_consequence=LOF_CSQ_SET,
            lof=lof_groups,
            mutation_type=["CpG"],
            transition=[True],
            cpg=[True],
            methylation_level=[0],
        ),
        "SFS_HCLoF_nonCpGtransition.txt.gz": hl.Struct(
            most_severe_consequence=LOF_CSQ_SET,
            lof=["HC"],
            mutation_type=["non-CpG transition"],
            transition=[True],
            cpg=[False],
            methylation_level=[0, 1, 2],
        ),
        "SFS_HCLoF_transversion.txt.gz": hl.Struct(
            most_severe_consequence=LOF_CSQ_SET,
            lof=["HC"],
            mutation_type=["transversion"],
            transition=[False],
            cpg=[False],
            methylation_level=[0, 1, 2],
        ),
        "SFS_HCLoF_methylCpG.txt.gz": hl.Struct(
            most_severe_consequence=LOF_CSQ_SET,
            lof=["HC"],
            mutation_type=["CpG"],
            transition=[True],
            cpg=[True],
            methylation_level=[2],
        ),
        "SFS_HCLoF_nonmethylCpG.txt.gz": hl.Struct(
            most_severe_consequence=LOF_CSQ_SET,
            lof=["HC"],
            mutation_type=["CpG"],
            transition=[True],
            cpg=[True],
            methylation_level=[0],
        ),
    }


def export_sfs_tables(
    agg_ht: hl.Table,
    output_dir: str,
    ms_csq_groups: Optional[List[str]] = None,
    group_by_ancestry: bool = False,
) -> None:
    """
    Export SFS tables for all combinations of functional categories and mutation types.

    :param agg_ht: Aggregated Hail table
    :param output_dir: Directory to save SFS tables
    :param ms_csq_groups: List of most severe consequence groups to process
    :param group_by_ancestry: Whether to export separate tables per ancestry group.
        Default is False.
    """
    if ms_csq_groups is None:
        ms_csq_groups = [
            "most_severe_consequence",
            "canonical_most_severe_consequence",
            "mane_select_most_severe_consequence",
        ]

    # Define AC groups.
    ac_groups = [
        "0",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "10-100",
        "100-1000",
        "1000-10000",
        ">10000",
    ]

    # Get grouping definitions.
    name_map = get_sfs_grouping_definitions()

    # Rename downsampling to Sample_size.
    agg_ht = agg_ht.annotate(**agg_ht.n).rename({"downsampling": "Sample_size"})

    # Check if gen_anc field exists in the table.
    has_gen_anc = "gen_anc" in agg_ht.row.dtype.fields

    if group_by_ancestry and not has_gen_anc:
        logger.warning(
            "group_by_ancestry is True but gen_anc field not found in aggregated table. "
            "Proceeding without ancestry grouping. "
            "Make sure to set --group-by-ancestry during aggregation step."
        )
        group_by_ancestry = False

    logger.info("Exporting SFS tables for %d combinations", len(name_map))

    for ms_csq in ms_csq_groups:
        logger.info("Processing most severe consequence group: %s", ms_csq)

        for file_name, grouping in name_map.items():
            logger.info("  Processing %s", file_name)

            if group_by_ancestry and has_gen_anc:
                # Get unique ancestry groups
                gen_ancs = sorted(
                    hl.eval(agg_ht.aggregate(hl.agg.collect_as_set(agg_ht.gen_anc)))
                )

                for gen_anc in gen_ancs:
                    logger.info("    Processing ancestry group: %s", gen_anc)

                    # Filter to this ancestry group
                    ht_anc = agg_ht.filter(agg_ht.gen_anc == gen_anc)

                    # Aggregate by sample size and AC group.
                    ht = ht_anc.group_by("Sample_size").aggregate(
                        _=hl.agg.group_by(
                            ht_anc.ac_group,
                            hl.agg.explode(
                                hl.agg.sum,
                                [
                                    ht_anc[ms_csq].get(
                                        hl.struct(
                                            most_severe_consequence=ms,
                                            lof=l,
                                            mutation_type=mt,
                                            transition=t,
                                            cpg=c,
                                            methylation_level=ml,
                                        )
                                    )
                                    for ms in grouping["most_severe_consequence"]
                                    for l in grouping["lof"]
                                    for mt in grouping["mutation_type"]
                                    for t in grouping["transition"]
                                    for c in grouping["cpg"]
                                    for ml in grouping["methylation_level"]
                                ],
                            ),
                        )
                    )

                    # Select AC groups as columns.
                    ht = ht.select(**{g: ht._.get(g, 0) for g in ac_groups})

                    # Export to file with ancestry group in path.
                    output_path = f"{output_dir}/{ms_csq}/{gen_anc}/{file_name}"
                    logger.info("      Exporting to %s", output_path)
                    ht.export(output_path)
            else:
                # Aggregate by sample size and AC group (aggregate across all
                # ancestries).
                ht = agg_ht.group_by("Sample_size").aggregate(
                    _=hl.agg.group_by(
                        agg_ht.ac_group,
                        hl.agg.explode(
                            hl.agg.sum,
                            [
                                agg_ht[ms_csq].get(
                                    hl.struct(
                                        most_severe_consequence=ms,
                                        lof=l,
                                        mutation_type=mt,
                                        transition=t,
                                        cpg=c,
                                        methylation_level=ml,
                                    )
                                )
                                for ms in grouping["most_severe_consequence"]
                                for l in grouping["lof"]
                                for mt in grouping["mutation_type"]
                                for t in grouping["transition"]
                                for c in grouping["cpg"]
                                for ml in grouping["methylation_level"]
                            ],
                        ),
                    )
                )

                # Select AC groups as columns.
                ht = ht.select(**{g: ht._.get(g, 0) for g in ac_groups})

                # Export to file.
                output_path = f"{output_dir}/{ms_csq}/{file_name}"
                logger.info("    Exporting to %s", output_path)
                ht.export(output_path)


def main(args):
    """Execute the SFS computation pipeline."""
    logger.info("Starting SFS computation pipeline")

    # Initialize Hail.
    hl.init(
        log="/sfs_computation.log",
        tmp_dir=args.tmp_dir,
    )

    genetic_ancestry_groups = (
        args.genetic_ancestry_groups if args.group_by_ancestry else None
    )

    # Step 1: Preprocess data.
    if args.preprocess_data:
        logger.info("Step 1: Preprocessing constraint data")
        preprocess_constraint_data(
            hl.read_table(args.input_ht_path), gen_ancs=genetic_ancestry_groups
        ).write(args.preprocessed_ht_path, overwrite=args.overwrite)

    # Step 2: Aggregate SFS data.
    if args.aggregate_data:
        logger.info("Step 2: Aggregating SFS data")
        aggregate_sfs_data(
            hl.read_table(args.preprocessed_ht_path),
            group_by_ancestry=args.group_by_ancestry,
        ).write(args.aggregated_ht_path, overwrite=args.overwrite)

    # Step 3: Export SFS tables.
    if args.export_tables:
        logger.info("Step 3: Exporting SFS tables")
        export_sfs_tables(
            hl.read_table(args.aggregated_ht_path),
            output_dir=args.output_dir,
            ms_csq_groups=args.ms_csq_groups,
            group_by_ancestry=args.group_by_ancestry,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute Site Frequency Spectrum (SFS) tables from constraint pipeline data"
    )

    # Input/Output paths
    parser.add_argument(
        "--input-ht-path",
        help="Path to input constraint table",
        type=str,
        default="gs://gnomad/v4.1/constraint_coverage_corrected/preprocessed_data/gnomad.v4.1.context.preprocessed.ht",
    )
    parser.add_argument(
        "--preprocessed-ht-path",
        help="Path to save preprocessed table",
        type=str,
        default="gs://gnomad-tmp-4day/julia/gnomad.v4.1.constraint.preprocessed.for_ac_group_computation.ht",
    )
    parser.add_argument(
        "--aggregated-ht-path",
        help="Path to save aggregated table",
        type=str,
        default="gs://gnomad-tmp-4day/julia/gnomad.v4.1.constraint.preprocessed.for_ac_group_computation.agg_ht.ht",
    )
    parser.add_argument(
        "--output-dir",
        help="Directory to save SFS tables",
        type=str,
        default="gs://gnomad-julia/loeuf_all/allele_count_by_downsampling",
    )
    parser.add_argument(
        "--tmp-dir",
        help="Temporary directory for Hail",
        type=str,
        default="gs://gnomad-tmp-4day",
    )

    # Pipeline steps
    parser.add_argument(
        "--preprocess-data",
        help="Preprocess constraint data",
        action="store_true",
    )
    parser.add_argument(
        "--aggregate-data",
        help="Aggregate SFS data",
        action="store_true",
    )
    parser.add_argument(
        "--export-tables",
        help="Export SFS tables",
        action="store_true",
    )

    # Options
    parser.add_argument(
        "--overwrite",
        help="Whether to overwrite existing output files",
        action="store_true",
    )
    parser.add_argument(
        "--ms-csq-groups",
        help="Most severe consequence groups to process",
        nargs="+",
        default=[
            "most_severe_consequence",
            "canonical_most_severe_consequence",
            "mane_select_most_severe_consequence",
        ],
    )
    parser.add_argument(
        "--genetic-ancestry-groups",
        help="Genetic ancestry groups to include in frequency data",
        nargs="+",
        choices=["afr", "amr", "eas", "nfe", "sas"],
        default=["afr", "amr", "eas", "nfe", "sas"],
    )
    parser.add_argument(
        "--group-by-ancestry",
        help="Whether to aggregate and export SFS tables separately by ancestry group",
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
