from typing import Union

import hail as hl

from .generic import *

def pre_process_data(ht: hl.Table, split_context_ht_path: str,
                     output_ht_path: str, overwrite: bool = False) -> None:
    """ 
    Add annotations from VEP context Table to gnomAD data.
    
    Function adds the following annotations:
        - context
        - methylation
        - coverage
        - gerp
        - pass_filters
    
    Function drops `a_index`, `was_split`, and`colocated_variants` annotations from gnomAD data.
    
    .. note::
        Function expects that multiallelic variants in the VEP context Table have been split.
    
    :param ht: gnomAD exomes or genomes public Hail Table. 
    :param split_context_ht_path: Path to VEP context Table.
    :param output_ht_path: Path to output Table.
    :param overwrite: Whether to overwrite existing data. Defaults to False.
    """
    context_ht = hl.read_table(split_context_ht_path).drop('a_index', 'was_split')
    context_ht = context_ht.annotate(vep=context_ht.vep.drop('colocated_variants'))
    ht.annotate(**context_ht[ht.key], pass_filters=hl.len(ht.filters) == 0).write(output_ht_path, overwrite)

def prepare_ht(ht, trimer: bool = False, annotate_coverage: bool = True):
    """
    Filter input Table and add annotations used in constraint calculations.
 
    Function filters to SNPs, removes rows with undefined contexts, collapses strands
    to deduplicate trimer or heptamer contexts, and annotates the input Table.
    
    The following annotations are added to the output Table:
        - ref
        - alt
        - methylation_level
        - exome_coverage
        - annotations added by `annotate_variant_types` and `collapse_strand`
   
    :param ht: Input Table to be annotated.
    :param trimer: Whether to use trimers or heptamers. Defaults to False.
    :param annotate_coverage: Whether to annotate the coverage of exome. Defaults to True.
    :return: Table with annotations.
    """
    str_len = 7
    
    if trimer:
        str_len = 3
        ht = trimer_from_heptamer(ht)

    if isinstance(ht, hl.Table):
        ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter((hl.len(ht.ref) == 1) & (hl.len(ht.alt) == 1) & ht.context.matches(f'[ATCG]{{{str_len}}}'))
        ht = annotate_variant_types(collapse_strand(ht), not trimer)
    else:
        ht = ht.annotate_rows(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter_rows((hl.len(ht.ref) == 1) & (hl.len(ht.alt) == 1) & ht.context.matches(f'[ATCG]{{{str_len}}}'))
        ht = annotate_variant_types(collapse_strand(ht), not trimer)
    annotation = {
        'methylation_level': hl.case().when(
            ht.cpg & (ht.methylation.MEAN > 0.6), 2
        ).when(
            ht.cpg & (ht.methylation.MEAN > 0.2), 1
        ).default(0)
    }
    if annotate_coverage:
        annotation['exome_coverage'] = ht.coverage.exomes.median
    return ht.annotate(**annotation) if isinstance(ht, hl.Table) else ht.annotate_rows(**annotation)
    
