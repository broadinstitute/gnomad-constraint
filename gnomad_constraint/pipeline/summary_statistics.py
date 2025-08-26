import argparse
import sys
from typing import List

import hail as hl
from gnomad.resources.grch38.gnomad import POPS
from gnomad.resources.grch38.reference_data import clinvar
from gnomad.utils.filtering import filter_low_conf_regions
from gnomad.utils.vep import (
    CSQ_ORDER,
    filter_vep_transcript_csqs,
    get_most_severe_consequence_for_summary,
)
from gnomad_qc.v3.resources.annotations import allele_data as v3_allele_data
from gnomad_qc.v4.resources.annotations import get_info

import gnomad_constraint.resources.resource_utils as constraint_res
from gnomad_constraint.utils.constraint import (
    annotate_mutation_type,
    count_variants_by_group,
)

EXOME_POPS = POPS["v4"]
GENOME_POPS = POPS["v3"]


root = "gs://gnomad-julia/v4.0/constraint"
subdir = "summary_results"
maps_ht_path = f"{root}/{subdir}/maps_plain_{{data_type}}.ht"
sfs_ht_path = f"{root}/{subdir}/sfs_{{data_type}}.ht"
loftee_maps_ht_path = f"{root}/{subdir}/maps_loftee_{{data_type}}.ht"
fifty_bp_maps_ht_path = f"{root}/{subdir}/maps_end_trunc_50bp_{{data_type}}.ht"
end_trunc_maps_ht_path = f"{root}/{subdir}/maps_end_trunc_gerp_{{data_type}}.ht"
loftee_assess_ht_path = f"{root}/{subdir}/freq_loftee_{{data_type}}.ht"
variants_per_sample_ht_path = f"{root}/{subdir}/variants_per_sample_{{data_type}}.ht"
observed_possible_ht_path = (
    f"{root}/{subdir}/observed_possible_expanded_{{data_type}}.ht"
)
observed_possible_sites_ht_path = (
    f"{root}/{subdir}/observed_possible_sites_{{data_type}}.txt"
)
indels_summary_ht_path = f"{root}/{subdir}/indels_summary_{{data_type}}.ht"
methylation_hist_file = f"{root}/{subdir}/methylation_hist.txt.bgz"
mutation_rate_ht_path = (
    f"{root}/{subdir}/mutation_rate_methylation_bins.v2_coverage_liftover.txt.bgz"
)
po_coverage_ht_path = f"{root}/{subdir}/prop_observed_by_coverage_no_common_pass_filtered_bins.v2_coverage_liftover.txt.bgz"
constraint_path_gene = (
    f"{root}/gnomad.v2.1.1.lof_metrics.by_gene.v2_coverage_liftover.txt.bgz"
)
constraint_path_transcript = (
    f"{root}/gnomad.v2.1.1.lof_metrics.by_transcript.v2_coverage_liftover.txt.bgz"
)


def get_downsamplings(ht):
    freq_meta = ht.freq_meta.collect()[0]
    downsamplings = [
        (i, int(x.get("downsampling")))
        for i, x in enumerate(freq_meta)
        if x.get("group") == "adj"
        and x.get("pop") == "global"
        and x.get("downsampling") is not None
    ]
    return downsamplings


def explode_downsamplings(ht, full_sample_size):
    downsamplings = get_downsamplings(ht)

    ht = ht.transmute(
        data=[
            hl.struct(
                downsampling=downsamplings[i][1],
                singletons=ht.singleton_downsampling_array[downsamplings[i][0]],
                observed=ht.downsampling_array[downsamplings[i][0]],
                possible=ht.possible,
            )
            for i in range(len(downsamplings))
        ]
        + [
            hl.struct(
                downsampling=full_sample_size,
                singletons=ht.singletons,
                observed=ht.observed,
                possible=ht.possible,
            )
        ]
    )
    ht = ht.explode("data")
    return ht.transmute(**ht.data)


def parse_lof_info(ht: hl.Table, location: hl.expr.StringExpression = None):
    """
    Must be exploded
    """
    if location is None:
        location = ht.vep.transcript_consequences.lof_info
    ht = ht.annotate(lof_data=location.split(",").map(lambda x: x.split(":")))
    return ht.transmute(
        lof_data=hl.dict(
            ht.lof_data.map(lambda x: (x[0], hl.if_else(hl.len(x) == 1, x[0], x[1])))
        )
    )


def get_bin_boundaries(
    ht: hl.Table, feature: str = "GERP_DIST", n_bins: int = 20, chrom: str = "10"
):
    temp_ht = hl.filter_intervals(ht, [hl.parse_locus_interval(chrom)])
    temp_ht = temp_ht.filter(hl.is_defined(temp_ht.lof_data.get(feature)))
    data = temp_ht.aggregate(hl.agg.collect(hl.float(temp_ht.lof_data.get(feature))))
    data = sorted(data)
    l = int(len(data) / n_bins)
    bin_boundaries = [data[x - 1] for x in range(len(data), -1, -l)]
    bin_boundaries.append(min(data))
    return bin_boundaries


def maps(
    ht: hl.Table,
    mutation_ht: hl.Table,
    additional_grouping: List[str] = [],
    singleton_expression: hl.expr.BooleanExpression = None,
    skip_worst_csq: bool = False,
) -> hl.Table:
    if not skip_worst_csq:
        additional_grouping.insert(0, "most_severe_csq")
    ht = count_variants_by_group(
        ht,
        count_singletons=True,
        additional_grouping=additional_grouping,
        use_table_group_by=True,
        singleton_expr=singleton_expression,
    )
    ht = ht.annotate(
        mu=mutation_ht[
            hl.struct(
                context=ht.context,
                ref=ht.ref,
                alt=ht.alt,
                methylation_level=ht.methylation_level,
            )
        ].mu_snp,
        ps=ht.singleton_count / ht.variant_count,
    )
    # if not ht.all(hl.is_defined(ht.mu)):
    #    print('Some mu were not found...')
    #    print(ht.aggregate(hl.agg.filter(hl.is_missing(ht.mu), hl.agg.take(ht.row, 1)[0])))
    #    sys.exit(1)
    syn_ps_ht = ht.filter(ht.most_severe_csq == "synonymous_variant")
    syn_ps_ht = syn_ps_ht.group_by(syn_ps_ht.mu).aggregate(
        singleton_count=hl.agg.sum(syn_ps_ht.singleton_count),
        variant_count=hl.agg.sum(syn_ps_ht.variant_count),
    )
    syn_ps_ht = syn_ps_ht.annotate(
        ps=syn_ps_ht.singleton_count / syn_ps_ht.variant_count
    )

    lm = syn_ps_ht.aggregate(
        hl.agg.linreg(
            syn_ps_ht.ps, [1, syn_ps_ht.mu], weight=syn_ps_ht.variant_count
        ).beta
    )
    print(f"Got MAPS calibration model of: slope: {lm[1]}, intercept: {lm[0]}")
    ht = ht.annotate(expected_singletons=(ht.mu * lm[1] + lm[0]) * ht.variant_count)

    agg_ht = ht.group_by(*additional_grouping).aggregate(
        singleton_count=hl.agg.sum(ht.singleton_count),
        expected_singletons=hl.agg.sum(ht.expected_singletons),
        variant_count=hl.agg.sum(ht.variant_count),
    )
    agg_ht = agg_ht.annotate(
        ps=agg_ht.singleton_count / agg_ht.variant_count,
        maps=(agg_ht.singleton_count - agg_ht.expected_singletons)
        / agg_ht.variant_count,
    )
    agg_ht = agg_ht.annotate(
        maps_sem=(agg_ht.ps * (1 - agg_ht.ps) / agg_ht.variant_count) ** 0.5
    )
    return agg_ht


def load_gtf_data():
    ht = hl.experimental.import_gtf(
        "gs://gnomad/resources/gencode/gencode.v39.annotation.gtf.gz",
        "GRCh38",
        True,
        min_partitions=12,
        force=True,
    )
    ht = ht.annotate(
        gene_id=ht.gene_id.split("\\.")[0],
        transcript_id=ht.transcript_id.split("\\.")[0],
        length=ht.interval.end.position - ht.interval.start.position + 1,
    )
    genes = (
        ht.filter(ht.feature == "gene")
        .select("gene_id", "gene_type", "gene_name", "length")
        .rename({"length": "gene_length"})
        .key_by("gene_id")
    )
    coding_regions = ht.filter(ht.feature == "CDS").select(
        "gene_id", "transcript_id", "transcript_type", "length", "level"
    )
    transcripts = (
        coding_regions.group_by(
            "transcript_id",
            "transcript_type",
            "gene_id",
            transcript_level=coding_regions.level,
        )
        .aggregate(
            cds_length=hl.agg.sum(coding_regions.length),
            num_coding_exons=hl.agg.count(),
        )
        .key_by("transcript_id")
    )
    return transcripts.annotate(**genes[transcripts.gene_id])


def select_primitives_from_ht(ht):
    return ht.select(
        **{
            x: v
            for x, v in ht.row_value.items()
            if v.dtype
            in {hl.tstr, hl.tint32, hl.tfloat32, hl.tint64, hl.tfloat64, hl.tbool}
        }
    )


def add_rank(ht, field, ascending=True, total_genes=None, bins=10, defined_only=False):
    if total_genes is None:
        if defined_only:
            total_genes = ht.aggregate(hl.agg.count_where(hl.is_defined(ht[field])))
        else:
            total_genes = ht.count()
    rank_field = ht[field] if ascending else -ht[field]
    ht = ht.order_by(rank_field).add_index(f"{field}_rank")
    ht = ht.annotate(
        **{
            f"{field}_rank": hl.or_missing(
                hl.is_defined(ht[field]), ht[f"{field}_rank"]
            )
        }
    )
    return ht.annotate(
        **{
            f"{field}_bin": hl.int(ht[f"{field}_rank"] * bins / total_genes),
            f"{field}_bin_6": hl.int(ht[f"{field}_rank"] * 6 / total_genes),
        }
    )


def main(args):
    version = "4.0"
    hl.init(
        log="/constraint_pipeline_summary_stats.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    if version.startswith("2"):
        build = "GRCh37"
    else:
        build = "GRCh38"

    if args.export_constraint_data:
        # ht = constraint_res.get_constraint_metrics_dataset(version=version).ht()
        ht = hl.read_table(
            "gs://gnomad/v4.0/constraint_liftover/metrics/gnomad.v4.0.constraint_metrics.ht"
        )
        gene_ht = load_gtf_data()
        gene = gene_ht.drop("gene_name")[ht.transcript]
        ht = ht.annotate(**gene)
        ht = ht.flatten()
        ht = ht.annotate(
            constraint_flag=hl.delimit(ht.constraint_flags, "|"),
            chromosome=ht.interval.start.contig,
            start_position=ht.interval.start.position,
            end_position=ht.interval.end.position,
        )
        ht = select_primitives_from_ht(ht)
        ht = ht.checkpoint(
            hl.utils.new_temp_file("constraint", "ht"), overwrite=args.overwrite
        )
        ht.export(constraint_path_transcript)

        # ht = hl.read_table("gs://gnomad-tmp-4day/constraint-Mw15brzVuepvLwZrBTRJCn.ht")
        ht = ht.filter(ht.mane_select)
        ht = (
            add_rank(ht, "lof.oe_ci.upper", defined_only=True)
            .drop("mane_select")
            .key_by("gene_id", "gene")
        )
        ht = select_primitives_from_ht(ht)
        ht.export(constraint_path_gene)

    if args.export_po_coverage:
        # po_ht = constraint_res.get_training_dataset(version, "autosome_par").ht()
        po_ht = hl.read_table(
            "gs://gnomad/v4.0/constraint_liftover/training_data/gnomad.v4.0.constraint_training.autosome_par.ht"
        )
        # po_ht = po_ht.filter(po_ht.exome_coverage <= 100)
        po_ht.describe()
        po_ht.show()
        po_ht.export(po_coverage_ht_path)

    # context_ht = constraint_res.vep_context_ht.versions[version].ht()
    # context_ht = constraint_res.get_preprocessed_ht("context", version, "autosome_par").ht()
    # context_ht = hl.filter_intervals(
    #    context_ht, [hl.parse_locus_interval("chrX", reference_genome=build)],
    #    keep=False
    # )
    mutation_ht = hl.read_table(
        "gs://gnomad/v4.0/constraint_liftover/mutation_rate/gnomad.v4.0.mutation_rate.ht"
    )  # constraint_res.get_mutation_ht(version).ht()

    data_type_sex_counts = {
        "exomes": {"XX": 367323, "XY": 363624},
        "genomes": {"XX": 38947, "XY": 37209},
    }
    if args.methylation_hist:
        context_ht = constraint_res.get_preprocessed_ht(
            "context", version, "autosome_par"
        ).ht()
        context_ht = hl.filter_intervals(
            context_ht,
            [hl.parse_locus_interval("chrX", reference_genome=build)],
            keep=False,
        )
        context_ht.describe()
        methylation_hist = context_ht.aggregate(
            hl.agg.hist(context_ht.methylation.methylation_level, 0, 15, 40)
        )
        data = list(zip(methylation_hist.bin_edges, methylation_hist.bin_freq))
        with hl.hadoop_open(methylation_hist_file, "w") as f:
            f.write("edge\tfreq\n")
            for edge, freq in data:
                f.write(f"{edge}\t{freq}\n")

    if args.export_mutation_rate:
        mutation_ht.export(mutation_rate_ht_path)

    for data_type, sample_sizes in data_type_sex_counts.items():
        print(f"Running {data_type}...")
        # ht = constraint_res.get_sites_resource(data_type, version).ht()
        ht = constraint_res.get_preprocessed_ht(data_type, version, "autosome_par").ht()
        ht = hl.filter_intervals(
            ht, [hl.parse_locus_interval("chrX", reference_genome=build)], keep=False
        )
        coverage_ht = constraint_res.get_coverage_ht(data_type, version).ht()

        # Obtain field name for median exome coverage.
        exome_median_cov_field = "median_approx"
        ht = filter_vep_transcript_csqs(
            ht,
            synonymous=False,
            filter_empty_csq=True,
            canonical=False,
            mane_select=True,
        )
        ht = get_most_severe_consequence_for_summary(ht)

        if args.run_indels:
            print(f"Running indels for {data_type}...")
            indel_ht = ht.filter(hl.is_indel(ht.alleles[0], ht.alleles[1]))
            indels = (
                indel_ht.group_by(
                    indel_ht.most_severe_csq,
                    indel_ht.coverage[data_type][exome_median_cov_field],
                    indel_length=hl.len(indel_ht.alleles[1])
                    - hl.len(indel_ht.alleles[0]),
                )
                .partition_hint(100)
                .aggregate(
                    singletons=hl.agg.count_where(indel_ht.freq[0].AC == 1),
                    observed=hl.agg.count(),
                    downsampling_array=hl.agg.array_sum(
                        indel_ht.freq.map(lambda x: x.AC > 0)
                    ),
                    singleton_downsampling_array=hl.agg.array_sum(
                        indel_ht.freq.map(lambda x: x.AC == 1)
                    ),
                )
            )
            total_real_estate = (
                coverage_ht.group_by(coverage_ht[exome_median_cov_field])
                .partition_hint(100)
                .aggregate(real_estate=hl.agg.count())
            )
            indels.annotate(
                coverage_real_estate=total_real_estate[
                    indels[exome_median_cov_field]
                ].real_estate
            ).write(indels_summary_ht_path.format(data_type=data_type), args.overwrite)
            indels = hl.read_table(indels_summary_ht_path.format(data_type=data_type))
            indels.export(
                indels_summary_ht_path.format(data_type=data_type).replace(
                    ".ht", ".txt.bgz"
                )
            )
            explode_downsamplings(
                indels.annotate(possible=indels.coverage_real_estate),
                sum(sample_sizes.values()),
            ).export(
                indels_summary_ht_path.format(data_type=data_type).replace(
                    ".ht", ".downsampling.txt.bgz"
                )
            )

        if args.run_sfs:
            criteria = hl.case(missing_false=True).when(ht.freq[0].AC == 1, "Singleton")
            if data_type == "genomes":
                criteria = criteria.when(ht.freq[0].AF < 1e-3, "< 0.1%")
            else:
                criteria = (
                    criteria.when(ht.freq[0].AC == 2, "Doubleton")
                    .when(ht.freq[0].AC <= 5, "AC 3-5")
                    .when(ht.freq[0].AF < 1e-4, "< 0.01%")
                    .when(ht.freq[0].AF < 1e-3, "0.01% - 0.1%")
                )
            sfs_ht = ht.annotate(
                freq_bin=criteria.when(ht.freq[0].AF < 1e-2, "0.1% - 1%")
                .when(ht.freq[0].AF < 1e-1, "1% - 10%")
                .default(">10%")
            )
            sfs_ht.group_by(
                sfs_ht.freq_bin,
                sfs_ht.most_severe_csq,
                snp=hl.is_snp(sfs_ht.alleles[0], sfs_ht.alleles[1]),
            ).aggregate(total=hl.agg.count()).write(
                sfs_ht_path.format(data_type=data_type), overwrite=args.overwrite
            )
            hl.read_table(sfs_ht_path.format(data_type=data_type)).export(
                sfs_ht_path.format(data_type=data_type).replace(".ht", ".txt.bgz")
            )
        snp_ht = ht
        if args.run_maps:
            print(f"Running MAPS for {data_type}...")
            maps_ht = maps(
                snp_ht,
                mutation_ht,
                additional_grouping=[
                    "protein_coding",
                ],
            )
            maps_ht.write(maps_ht_path.format(data_type=data_type), args.overwrite)
            hl.read_table(maps_ht_path.format(data_type=data_type)).export(
                maps_ht_path.format(data_type=data_type).replace(".ht", ".txt.bgz")
            )

        if args.run_loftee_maps:
            exploded_snp_ht = snp_ht.transmute(
                transcript_consequences=snp_ht.vep.transcript_consequences
            )
            exploded_snp_ht = exploded_snp_ht.explode(
                exploded_snp_ht.transcript_consequences
            )
            exploded_snp_ht = parse_lof_info(
                exploded_snp_ht, exploded_snp_ht.transcript_consequences.lof_info
            )
            exploded_snp_ht = exploded_snp_ht.annotate(
                most_severe_csq=hl.literal(CSQ_ORDER).find(
                    lambda x: exploded_snp_ht.transcript_consequences.consequence_terms.contains(
                        x
                    )
                ),
                lof=exploded_snp_ht.transcript_consequences.lof,
                lof_filter=exploded_snp_ht.transcript_consequences.lof_filter,
                lof_flags=exploded_snp_ht.transcript_consequences.lof_flags,
                lof_filter_simplified=hl.if_else(
                    exploded_snp_ht.transcript_consequences.lof_filter.contains(","),
                    "MULTIPLE",
                    exploded_snp_ht.transcript_consequences.lof_filter,
                ),
            )

            maps_ht = maps(
                exploded_snp_ht, mutation_ht, ["lof", "lof_filter", "lof_flags"]
            )
            maps_ht.write(
                loftee_maps_ht_path.format(data_type=data_type), args.overwrite
            )
            hl.read_table(loftee_maps_ht_path.format(data_type=data_type)).export(
                loftee_maps_ht_path.format(data_type=data_type).replace(
                    ".ht", ".txt.bgz"
                )
            )

        if args.run_end_trunc_maps:
            exploded_snp_ht = snp_ht.transmute(
                transcript_consequences=snp_ht.vep.transcript_consequences
            )
            exploded_snp_ht = exploded_snp_ht.explode(
                exploded_snp_ht.transcript_consequences
            )
            exploded_snp_ht = parse_lof_info(
                exploded_snp_ht, exploded_snp_ht.transcript_consequences.lof_info
            )

            gerp_bin_boundaries = get_bin_boundaries(
                exploded_snp_ht, "GERP_DIST", n_bins=40
            )
            bp_bin_boundaries = get_bin_boundaries(
                exploded_snp_ht, "BP_DIST", n_bins=40
            )
            percentile_breaks = list(map(lambda x: x / 10, range(9, 0, -1))) + list(
                map(lambda x: x / 200, range(18, -1, -1))
            )

            exploded_snp_ht = exploded_snp_ht.annotate(
                most_severe_csq=hl.literal(CSQ_ORDER).find(
                    lambda x: exploded_snp_ht.transcript_consequences.consequence_terms.contains(
                        x
                    )
                ),
                fifty_bp_rule=exploded_snp_ht.lof_data.get("50_BP_RULE"),
                gerp=hl.literal(gerp_bin_boundaries).find(
                    lambda x: hl.float(exploded_snp_ht.lof_data.get("GERP_DIST")) >= x
                ),
                bp=hl.literal(bp_bin_boundaries).find(
                    lambda x: hl.float(exploded_snp_ht.lof_data.get("BP_DIST")) >= x
                ),
                percentile=hl.literal(percentile_breaks).find(
                    lambda x: hl.float(exploded_snp_ht.lof_data.get("PERCENTILE")) >= x
                ),
            )
            # maps_ht = maps(exploded_snp_ht, mutation_ht, ['bp', ])
            # maps_ht = maps(exploded_snp_ht, mutation_ht, ['percentile', ])

            maps_ht = maps(
                exploded_snp_ht,
                mutation_ht,
                [
                    "fifty_bp_rule",
                ],
            )
            maps_ht.write(
                fifty_bp_maps_ht_path.format(data_type=data_type), args.overwrite
            )
            hl.read_table(fifty_bp_maps_ht_path.format(data_type=data_type)).export(
                fifty_bp_maps_ht_path.format(data_type=data_type).replace(
                    ".ht", ".txt.bgz"
                )
            )
            exploded_snp_ht = exploded_snp_ht.filter(
                (exploded_snp_ht.most_severe_csq == "synonymous_variant")
                | (exploded_snp_ht.fifty_bp_rule == "FAIL")
            )
            maps_ht = maps(
                exploded_snp_ht,
                mutation_ht,
                [
                    "gerp",
                ],
            )
            maps_ht.write(
                end_trunc_maps_ht_path.format(data_type=data_type), args.overwrite
            )
            hl.read_table(end_trunc_maps_ht_path.format(data_type=data_type)).export(
                end_trunc_maps_ht_path.format(data_type=data_type).replace(
                    ".ht", ".txt.bgz"
                )
            )

        if args.run_obs_poss:
            print(f"Running observed/possible for {data_type}...")

            context_ht = constraint_res.get_preprocessed_ht(
                "context", version, "autosome_par"
            )
            context_ht = filter_vep_transcript_csqs(
                context_ht,
                synonymous=False,
                filter_empty_csq=True,
                canonical=False,
                mane_select=True,
            )
            context_ht = get_most_severe_consequence_for_summary(context_ht)

            observed = (
                snp_ht.group_by(
                    snp_ht.most_severe_csq,
                    snp_ht.context,
                    snp_ht.ref,
                    snp_ht.alt,
                    snp_ht.methylation_level,
                    snp_ht.coverage[data_type][exome_median_cov_field],
                )
                .partition_hint(1000)
                .aggregate(
                    observed=hl.agg.count(),
                    singletons=hl.agg.count_where(snp_ht.freq[0].AC == 1),
                    downsampling_array=hl.agg.array_sum(
                        snp_ht.freq.map(lambda x: x.AC > 0)
                    ),
                    singleton_downsampling_array=hl.agg.array_sum(
                        snp_ht.freq.map(lambda x: x.AC == 1)
                    ),
                )
            )
            possible = (
                context_ht.group_by(
                    context_ht.most_severe_csq,
                    context_ht.context,
                    context_ht.ref,
                    context_ht.alt,
                    context_ht.methylation_level,
                    context_ht.coverage[data_type][exome_median_cov_field],
                )
                .partition_hint(1000)
                .aggregate(possible=hl.agg.count())
            )
            annotate_mutation_type(observed.join(possible), False).write(
                observed_possible_ht_path.format(data_type=data_type), args.overwrite
            )

            explode_downsamplings(
                hl.read_table(observed_possible_ht_path.format(data_type=data_type)),
                sum(sample_sizes.values()),
            ).export(
                observed_possible_ht_path.format(data_type=data_type).replace(
                    ".ht", ".txt.bgz"
                )
            )

        if args.run_obs_poss_sites:
            print(f"Running observed/possible sites for {data_type}...")
            if data_type == "exomes":
                initial_alleles = get_info().ht().select("allele_data")
            else:
                initial_alleles = v3_allele_data.ht()
            initial_alleles = initial_alleles.filter(initial_alleles.a_index == 1)

            sites_snp_ht = snp_ht.filter(hl.is_defined(initial_alleles[snp_ht.key]))

            observed = sites_snp_ht.aggregate(
                hl.agg.counter(sites_snp_ht.coverage[data_type][exome_median_cov_field])
            )
            possible = context_ht.aggregate(
                hl.agg.counter(
                    context_ht.coverage[data_type][data_type][exome_median_cov_field]
                )
            )
            for coverage, n_sites in observed.items():
                if coverage not in possible:
                    raise ValueError(f"{coverage} not found in possible dict")

            temp_f = hl.utils.new_local_temp_file()
            print(f"Writing to {temp_f}...")
            with open(temp_f, "w") as f:
                f.write("coverage\tn_sites_observed\tn_sites_possible\n")
                for coverage, n_sites_possible in possible.items():
                    n_sites_observed = observed.get(coverage, 0)
                    f.write(f"{coverage}\t{n_sites_observed}\t{n_sites_possible / 3}\n")
            hl.hadoop_copy(
                f"file://{temp_f}",
                observed_possible_sites_ht_path.format(data_type=data_type),
            )

        if data_type == "exomes" and args.assess_loftee:
            ht = ht.annotate(
                vep=ht.vep.annotate(
                    transcript_consequences=ht.vep.transcript_consequences.filter(
                        lambda x: hl.is_defined(x.lof)
                    )
                )
            )

            ht = ht.annotate(
                freq_bin=hl.case(missing_false=True)
                .when(ht.freq[0].AC == 1, "Singleton")
                .when(ht.freq[0].AC == 2, "Doubleton")
                .when(ht.freq[0].AF < 1e-4, "< 0.01%")
                .when(ht.freq[0].AF < 1e-3, "0.01% - 0.1%")
                .when(ht.freq[0].AF < 1e-2, "0.1% - 1%")
                .when(ht.freq[0].AF < 1e-1, "1% - 10%")
                .default(">10%"),
                fail_loftee=hl.any(
                    lambda tc: (tc.lof == "LC"), ht.vep.transcript_consequences
                ),
                loftee_os=hl.any(
                    lambda tc: (tc.lof == "OS"), ht.vep.transcript_consequences
                ),
                pass_loftee=hl.any(
                    lambda tc: (tc.lof == "HC"), ht.vep.transcript_consequences
                ),
                pass_loftee_with_flags=hl.any(
                    lambda tc: (tc.lof == "HC") & hl.is_missing(tc.lof_flags),
                    ht.vep.transcript_consequences,
                ),
            )
            ht = ht.group_by(
                ht.freq_bin,
                ht.most_severe_csq,
                ht.fail_loftee,
                ht.loftee_os,
                ht.pass_loftee,
                ht.pass_loftee_with_flags,
            ).aggregate(n=hl.agg.count())

            clinvar_ht = clinvar.ht()
            path_labels = hl.literal(
                {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
            )
            clinvar_ht = clinvar_ht.filter(
                hl.or_missing(
                    hl.len(clinvar_ht.info.CLNSIG) > 0,
                    path_labels.contains(clinvar_ht.info.CLNSIG[0]),
                )
                & hl.is_missing(clinvar_ht.info.CLNSIGCONF)
            )
            clinvar_ht = filter_vep_transcript_csqs(
                clinvar_ht,
                synonymous=False,
                filter_empty_csq=True,
                canonical=False,
                mane_select=True,
            )
            clinvar_ht = get_most_severe_consequence_for_summary(clinvar_ht)

            clinvar_ht = clinvar_ht.annotate(
                vep=clinvar_ht.vep.annotate(
                    transcript_consequences=clinvar_ht.vep.transcript_consequences.filter(
                        lambda x: hl.is_defined(x.lof)
                    )
                )
            )
            clinvar_ht = clinvar_ht.annotate(
                freq_bin=hl.delimit(clinvar_ht.info.CLNREVSTAT[:2], ",").replace(
                    "_", " "
                ),
                fail_loftee=hl.any(
                    lambda tc: (tc.lof == "LC"), clinvar_ht.vep.transcript_consequences
                ),
                loftee_os=hl.any(
                    lambda tc: (tc.lof == "OS"), clinvar_ht.vep.transcript_consequences
                ),
                pass_loftee=hl.any(
                    lambda tc: (tc.lof == "HC"), clinvar_ht.vep.transcript_consequences
                ),
                pass_loftee_with_flags=hl.any(
                    lambda tc: (tc.lof == "HC") & hl.is_missing(tc.lof_flags),
                    clinvar_ht.vep.transcript_consequences,
                ),
            )
            clinvar_ht = clinvar_ht.group_by(
                clinvar_ht.freq_bin,
                clinvar_ht.most_severe_csq,
                clinvar_ht.fail_loftee,
                clinvar_ht.loftee_os,
                clinvar_ht.pass_loftee,
                clinvar_ht.pass_loftee_with_flags,
            ).aggregate(n=hl.agg.count())

            ht = ht.union(clinvar_ht)

            ht.write(loftee_assess_ht_path.format(data_type=data_type), args.overwrite)
            hl.read_table(loftee_assess_ht_path.format(data_type=data_type)).export(
                loftee_assess_ht_path.format(data_type=data_type).replace(
                    ".ht", ".txt.bgz"
                )
            )

        if args.run_variants_per_sample:
            ht = constraint_res.get_sites_resource(data_type, version).ht()

            ht = ht.filter(
                (hl.len(ht.filters) == 0)
            )  # & get_an_adj_criteria(ht, sample_sizes))
            ht = filter_vep_transcript_csqs(
                ht,
                synonymous=False,
                filter_empty_csq=True,
                canonical=False,
                mane_select=True,
            )
            ht = get_most_severe_consequence_for_summary(ht)
            ht = filter_low_conf_regions(ht, annotate_instead_of_filter=True)

            def build_criteria(ht: hl.Table, data_type: str, index: int = 0):
                criteria = (
                    hl.case(missing_false=True)
                    .when(ht.freq[index].AC == 0, "Not found")
                    .when(ht.freq[index].AC == 1, "Singleton")
                    .when(ht.freq[index].AC == 2, "Doubleton")
                )
                if data_type == "genomes":
                    criteria = criteria.when(ht.freq[index].AF < 1e-3, "< 0.1%")
                else:
                    criteria = (
                        criteria.when(ht.freq[index].AC <= 5, "AC 3 - 5")
                        .when(ht.freq[index].AF < 1e-4, "AC 6 - 0.01%")
                        .when(ht.freq[index].AF < 1e-3, "0.01% - 0.1%")
                    )
                return (
                    criteria.when(ht.freq[index].AF < 1e-2, "0.1% - 1%")
                    .when(ht.freq[index].AF < 1e-1, "1% - 10%")
                    .when(ht.freq[index].AF > 0.95, ">95%")
                    .default("10% - 95%")
                )

            pops = list(
                map(
                    lambda x: x.lower(),
                    EXOME_POPS if data_type == "exomes" else GENOME_POPS,
                )
            )
            pops = [
                (pop, hl.eval(ht.freq_index_dict[f"gnomad_{pop}"])) for pop in pops
            ] + [("global", 0)]
            ht = ht.group_by(
                "most_severe_csq",
                "lof",
                "no_lof_flags",
                "protein_coding",
                "curated",
                **ht.regions,
            ).aggregate(
                pop_bin_sums=[
                    (
                        pop,
                        hl.agg.group_by(
                            build_criteria(ht, data_type, index),
                            [
                                hl.agg.sum(ht.freq[index].AC),
                                hl.agg.sum(ht.freq[index].homozygote_count),
                            ],
                        ),
                    )
                    for pop, index in pops
                ]
            )
            ht = ht.explode("pop_bin_sums")
            ht = ht.transmute(
                pop=ht.pop_bin_sums[0], bin_sums=hl.array(ht.pop_bin_sums[1])
            )
            ht = ht.explode("bin_sums")
            ht = ht.transmute(
                bin=ht.bin_sums[0], total=ht.bin_sums[1][0], total_hom=ht.bin_sums[1][1]
            )
            # # Alternative to above code:
            # pop_freq_mapping = {f'bin_{pop}': build_criteria(ht, data_type, index) for pop, index in pops}
            # ht = ht.annotate(**pop_freq_mapping)
            # ht = ht.group_by(*list(pop_freq_mapping), 'worst_csq', 'lof', 'no_lof_flags', 'protein_coding').aggregate(
            #     **{f'count_{pop}': hl.agg.sum(ht.freq[index].AC) for pop, index in pops},
            #     ** {f'hom_{pop}': hl.agg.sum(ht.freq[index].homozygote_count) for pop, index in pops}
            # )
            ht.write(
                variants_per_sample_ht_path.format(data_type=data_type),
                overwrite=args.overwrite,
            )
            hl.read_table(
                variants_per_sample_ht_path.format(data_type=data_type)
            ).export(
                variants_per_sample_ht_path.format(data_type=data_type).replace(
                    ".ht", ".txt.bgz"
                )
            )

    if args.print_model_info:
        ht = constraint_res.get_predicted_proportion_observed_dataset(
            version=version
        ).ht()
        print(hl.eval(ht.apply_model_params))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--overwrite", help="Overwrite everything", action="store_true")
    parser.add_argument(
        "--run_indels", help="Overwrite everything", action="store_true"
    )
    parser.add_argument("--run_sfs", help="Overwrite everything", action="store_true")
    parser.add_argument("--run_maps", help="Overwrite everything", action="store_true")
    parser.add_argument(
        "--run_loftee_maps", help="Overwrite everything", action="store_true"
    )
    parser.add_argument(
        "--run_end_trunc_maps", help="Overwrite everything", action="store_true"
    )
    parser.add_argument(
        "--run_obs_poss", help="Overwrite everything", action="store_true"
    )
    parser.add_argument(
        "--run_obs_poss_sites", help="Overwrite everything", action="store_true"
    )
    parser.add_argument(
        "--run_variants_per_sample", help="Overwrite everything", action="store_true"
    )
    parser.add_argument(
        "--assess_loftee", help="Overwrite everything", action="store_true"
    )
    parser.add_argument(
        "--methylation_hist", help="Overwrite everything", action="store_true"
    )
    parser.add_argument(
        "--export_mutation_rate", help="Overwrite everything", action="store_true"
    )
    parser.add_argument(
        "--export_po_coverage", help="Overwrite everything", action="store_true"
    )
    parser.add_argument(
        "--print_model_info", help="Overwrite everything", action="store_true"
    )
    parser.add_argument(
        "--export_constraint_data", help="Overwrite everything", action="store_true"
    )

    args = parser.parse_args()

    main(args)
