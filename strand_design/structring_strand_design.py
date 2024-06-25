import ipy_oxdna as iox
import ipy_oxdna.dna_structure as dna
import numpy as np
import argparse
from pathlib import Path

def parse_arguments():
    parser = argparse.ArgumentParser(description='FATCAT and TM-Align Strucutral Homology Screens')
    parser.add_argument('-s', '--strucutre_file', metavar='strucutre_file', type=str, help='Directory containing the query PDB files EX: ./query_dir/')
    parser.add_argument('-i', '--input_md_file', metavar='input_md_file', type=str, help='Directory of proteome directories EX: ./proteomes/ (where ./proteomes contains multiple dirs called i.e human, mouse, drome)')
    parser.add_argument('-c', '--coding_sequence', metavar='output_name', type=str, help='Location to save the csv formatted results of the FATCAT and TM-Align search. Default: ./fatcat_tmalign_homology_search.csv')
    parser.add_argument('-o', '--output_name', metavar='output_name', type=str, help='Location to save the csv formatted results of the FATCAT and TM-Align search. Default: ./fatcat_tmalign_homology_search.csv')
    parser.add_argument('--force_overwrite', action='store_true', help='Force overwrite of existing results if they exist.')
    args = parser.parse_args()
    return args, parser

def main():
    strucutre_file, input_md_file, coding_sequence, output_name = collect_args()
    


def collect_args():
    """
    Checks if the provided arguments exist and returns the path to:
    strucutre_file, input_md_file, coding_sequence, and output_name

    Raises:
        ValueError: Raises value error if the provided arguments are not valid

    Returns:
        The path to the strucutre_file, input_md_file, coding_sequence, and output_name
    """
    args, parser = parse_arguments()

    if not args.strucutre_file:
        parser.print_help()
        raise ValueError('Please provide a strucutre_file using the -s or --strucutre_file flag')
    else:
        strucutre_file = Path(args.strucutre_file)
        if not strucutre_file.exists():
            raise ValueError(f'Provided strucutre_file path {strucutre_file} does not exist')

    if not args.input_md_file:
        parser.print_help()
        raise ValueError('Please provide a input_md_file using the -i or --input_md_file flag')
    else:
        input_md_file = Path(args.input_md_file)
        if not input_md_file.exists():
            raise ValueError(f'Provided input_md_file path {input_md_file} does not exist')

    if not args.coding_sequence:
        parser.print_help()
        raise ValueError('Please provide a coding_sequence using the -c or --coding_sequence flag')
    else:
        coding_sequence = Path(args.coding_sequence)
        if not coding_sequence.exists():
            raise ValueError(f'Provided coding_sequence path {coding_sequence} does not exist')

    if not args.output_name:
        output_name = 'coding_sequence_embbeded_strucutre'
    else:
        output_name = args.output_name

    if not args.force_overwrite:
        if (Path(f'{output_name}.dat').exists()) or (Path(f'{output_name}.top').exists()):
            raise ValueError(f'Output file {output_name} already exists.\
                Please use the --force_overwrite flag to overwrite the file.')
    return strucutre_file, input_md_file, coding_sequence, output_name

def get_structure():
    pass


def get_hblist():
    pass


def get_halfs():
    pass


def get_idx_of_coding_complement():
    pass


def find_where_coding_binds_coding():
    pass


def mutate_strucutre():
    pass


def mutate_coding():
    pass


def mutate_noncoding():
    pass


def mutate_leftover_non_mutated_coding():
    pass


if __name__ == '__main__':
    main()