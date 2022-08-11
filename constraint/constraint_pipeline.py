import hail as hl

from gnomad.resources.grch38.reference_data import vep_context 
from gnomad_qc.v2.resources.basics import(
    get_gnomad_public_data
)
from utils import *

context_ht_path = vep_context.versions['101'].path

def main(args):
    if args.pre_process_data:
        pre_process_data(get_gnomad_public_data('genomes'), context_ht_path, processed_genomes_ht_path, args.overwrite)
        pre_process_data(get_gnomad_public_data('exomes'), context_ht_path, processed_exomes_ht_path, args.overwrite)

    full_context_ht = prepare_ht(hl.read_table(context_ht_path), args.trimers)
    full_genome_ht = prepare_ht(hl.read_table(processed_genomes_ht_path), args.trimers)
    full_exome_ht = prepare_ht(hl.read_table(processed_exomes_ht_path), args.trimers)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--trimers', help='Use trimers instead of heptamers', action='store_true')
    parser.add_argument('--pre_process_data', help='Pre-process data', action='store_true')
