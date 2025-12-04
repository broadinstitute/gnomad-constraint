import hail as hl
import argparse
import logging

from gnomad_constraint.experimental.proemis3d.genetics_gym.constants import (
    KEY_GROUPS,
    UNIPROT_ID_FIELD,
    LINKER_KEY_FIELDS,
    SCORE_KEY_GROUPS,
    SCORE_FIELDS,
    HIGHER_IS_LESS_DELETERIOUS,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("genetics_gym_missense_scores_combine")
logger.setLevel(logging.INFO)

BASE_HT_PATH = "gs://gnomad-julia/genetics_gym/linkers/vsm_base_em_gene.ht"
KEYED_HT_PATH = "gs://gnomad-tmp-4day/genetics_gym/linkers/vsm_base_em_gene_keyed.ht"
PARTITIONS_HE_PATH = "gs://gnomad-julia/genetics_gym/linkers/vsm_base_em_gene_partitions.he"


def main(args):
    """Combine missense scores for Genetics Gym."""
    hl.init(
        log="/genetics_gym_missense_scores_combine.log",
        #tmp_dir="gs://trisha-tmp",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    overwrite = args.overwrite
    test = args.test

    if args.rekey_and_partition_base_ht:
        logger.info("Rekeying and partitioning base HT...")
        ht = hl.read_table(BASE_HT_PATH)
        ht = ht.transmute(**{UNIPROT_ID_FIELD: ht.uniprot_em_id})
        tmp_path = hl.utils.new_temp_file("vsm_base_em_gene_keyed", "ht")
        ht = ht.key_by(*LINKER_KEY_FIELDS).checkpoint(tmp_path)
        hl.experimental.write_expression(
            ht._calculate_new_partitions(1000), PARTITIONS_HE_PATH, overwrite=overwrite
        )
        partition_intervals = hl.experimental.read_expression(PARTITIONS_HE_PATH)
        ht = hl.read_table(tmp_path, _intervals=hl.eval(partition_intervals))
        ht = ht.checkpoint(KEYED_HT_PATH, overwrite=overwrite)
        logger.info("Finished rekeying and partitioning base HT.")
    
    if args.create_subset_key_partitions:
        logger.info("Creating subset key partitions...")
        partition_intervals = hl.experimental.read_expression(PARTITIONS_HE_PATH)
        for key_name, keys in KEY_GROUPS.items():
            key_group_partitions = [
                hl.utils.interval.Interval(
                    start=interval.start.select(*keys),
                    end=interval.end.select(*keys),
                    includes_start=True, 
                    includes_end=False,
                    point_type=partition_intervals[0].start.select(*keys).dtype,
                )
                for interval in hl.eval(partition_intervals)
            ]
            hl.experimental.write_expression(
                key_group_partitions, 
                f"gs://gnomad-julia/genetics_gym/linkers/vsm_base_{key_name}_partitions.he", 
                overwrite=overwrite
            )
        logger.info("Finished creating subset key partitions.")
    
    if args.combine_scores:
        logger.info("Combining scores...")
        base_ht = hl.read_table(KEYED_HT_PATH)
        if test:
            base_ht = base_ht._filter_partitions(range(5))

        scores_expr = {}
        for score, key_group in SCORE_KEY_GROUPS.items():
            partition_intervals = hl.experimental.read_expression(
                f"gs://gnomad-julia/genetics_gym/linkers/vsm_base_{key_group}_partitions.he"
            )
            ht = hl.read_table(
                f"gs://gnomad-julia/genetics_gym/imported_vsm/{score}.ht",
                _intervals=hl.eval(partition_intervals)
            )

            if test:
                ht = ht._filter_partitions(range(5))

            ht_indexed = ht.index(*[base_ht[k] for k in KEY_GROUPS[key_group]])

            for s, add_neg in HIGHER_IS_LESS_DELETERIOUS[score].items():
                scores_expr[s] = ht_indexed[s]
                if add_neg:
                    scores_expr[f"{s}_neg"] = -ht_indexed[s]
        
        base_ht = base_ht.annotate(**scores_expr)

        base_ht.describe()
        base_ht = base_ht.checkpoint(
            "gs://gnomad-tmp-4day/genetics_gym/linkers/vsm_base_em_gene.annotated_scores.ht",
            overwrite=overwrite
        )
        base_ht.show()

        logger.info("Finished combining scores.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--overwrite", help="Whether to overwrite output files.", action="store_true"
    )
    parser.add_argument(
        "--test",
        help="Whether to run a test instead of the full pipeline",
        action="store_true",
    )
    parser.add_argument(
        "--rekey-and-partition-base-ht", 
        help="Rekey and partition the base HT.",
        action="store_true"
    )
    parser.add_argument(
        "--create-subset-key-partitions",
        help="Create subset key partitions.",
        action="store_true"
    )
    parser.add_argument(
        "--combine-scores",
        help="Combine scores.",
        action="store_true"
    )

    args = parser.parse_args()
    main(args)
