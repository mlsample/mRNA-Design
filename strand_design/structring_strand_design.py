import ipy_oxdna as iox
import ipy_oxdna.dna_structure as dna
from oxDNA_analysis_tools.output_bonds import output_bonds
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, inbox, describe
from oxDNA_analysis_tools.UTILS.oxview import oxdna_conf
import numpy as np
import argparse
import subprocess
from pathlib import Path
import sys
import pandas as pd

def main():
    strucutre_file, topology_file, input_md_file, coding_sequence, output_name = collect_args()

    coding_seq_str = parse_coding_sequence(coding_sequence)

    strucutre = get_structure(strucutre_file, topology_file)

    hb_id_1, hb_id_2 = get_hblist(strucutre_file, topology_file, input_md_file)
    
    half, indexes_1, indexes_2 = get_halfs(strucutre)
    
    pair_map = make_pair_map(hb_id_1, hb_id_2)
    index_to_seq_map = get_index_to_seq_map(strucutre)
    
    coding_indexes, coding_complement = get_idx_of_coding_complement(indexes_1, indexes_2, pair_map, coding_seq_str)
    
    coding_binds_coding_indexes = find_where_coding_binds_coding(coding_indexes, coding_complement, pair_map)
    
    return 0


def parse_arguments():
    parser = argparse.ArgumentParser(description='Embeds a coding sequence into a ssOrigami strucutre file')
    parser.add_argument('-s', '--strucutre_file', metavar='strucutre_file', type=str,
                        help='Path to oxDNA .dat file of ssOrigami strucutre EX: ./ssOrigami.dat')
    parser.add_argument('-t', '--topology_file', metavar='topology_file', type=str,
                        help='Path to oxDNA .top file of ssOrigami file EX: ./ssOrigami.dat')
    parser.add_argument('-i', '--input_md_file', metavar='input_md_file', type=str,
                        help='Path to oxDNA input file EX: ./input')
    parser.add_argument('-c', '--coding_sequence', metavar='output_name', type=str,
                        help='Path to .txt or .fasta file coding sequence EX: ./coding_seq.txt')
    parser.add_argument('-o', '--output_name', metavar='output_name', type=str,
                        help='name of the Default: ./fatcat_tmalign_homology_search.csv')
    parser.add_argument('--force_overwrite', action='store_true',
                        help='Force overwrite of existing results if they exist.')

    args = parser.parse_args()
    return args, parser

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
        assert strucutre_file.suffix == '.dat', f"strucutre_file must be a .dat file and it currently {strucutre_file.suffix}"


    if not args.topology_file:
        parser.print_help()
        raise ValueError('Please provide a topology_file using the -t or --topology_file flag')
    else:
        topology_file = Path(args.topology_file)
        if not topology_file.exists():
            raise ValueError(f'Provided topology_file path {topology_file} does not exist')
        assert topology_file.suffix == '.top', "topology_file must be a .top file"

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
        if coding_sequence.suffix not in ['.txt', '.fasta']:
            raise ValueError('coding_sequence must be a .txt or .fasta file')
        with open(coding_sequence, 'r') as f:
            lines = f.readlines()
            if coding_sequence.suffix == '.fasta':
                assert len(lines) == 2, "fasta file must have 2 lines"
            else:
                assert len(lines) == 1, "txt file must have 1 line"

    if not args.output_name:
        output_name = 'coding_sequence_embbeded_strucutre'
    else:
        output_name = args.output_name

    if not args.force_overwrite:
        if (Path(f'{output_name}.dat').exists()) or (Path(f'{output_name}.top').exists()):
            raise ValueError(f'Output file {output_name} already exists.\
                Please use the --force_overwrite flag to overwrite the file.')
    return strucutre_file, topology_file, input_md_file, coding_sequence, output_name


def parse_coding_sequence(coding_sequence):
    if coding_sequence.suffix == '.fasta':
        with open(coding_sequence, 'r') as f:
            lines = f.readlines()
            coding_seq_str = lines[1]
    elif coding_sequence.suffix == '.txt':
        with open(coding_sequence, 'r') as f:
            coding_seq_str = f.readline()[0]

    unique_bases = set(coding_seq_str)
    possible_bases = set('ATCGUatcgu')

    assert unique_bases.issubset(possible_bases), "coding_sequence must only contain A, T, C, G, or U"
    coding_seq_str = coding_seq_str.upper()
    return coding_seq_str


def get_structure(strucutre_file, topology_file):
    strucutre = dna.load_dna_structure(topology_file, strucutre_file)
    return strucutre


def get_hblist(strucutre_file, topology_file, input_md_file):

    inline_code =  f"""
from oxDNA_analysis_tools.output_bonds import output_bonds
from oxDNA_analysis_tools.UTILS.RyeReader import describe
top_info, traj_info = describe("{topology_file.as_posix()}", "{strucutre_file.as_posix()}")
output_bonds(traj_info, top_info, "{input_md_file.as_posix()}")
    """

    result = subprocess.run([sys.executable, '-c', inline_code], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)

    stdout_output = result.stdout

    potential_information = [info for info in stdout_output.split('\n') if "#" not in info]
    parsed_potentials = [info.strip().split(' ') for info in potential_information]
    parsed_potentials = np.array([info for info in parsed_potentials if len(info) == 11], dtype=float)

    headers = stdout_output.split('\n')[0].split(' ')
    potential_names = headers[:-3]
    potential_names[0] = potential_names[0][1:]
    potential_names[-1] = potential_names[-1][:-1]
    
    df = pd.DataFrame(parsed_potentials, columns=potential_names)
    hb_series = df[df['HB'] < -0.1]
    hb_id_1 = hb_series['id1']
    hb_id_2 = hb_series['id2']
    
    return hb_id_1, hb_id_2


def get_halfs(strucutre):
    n_bases = strucutre.get_num_bases()
    half = int(np.floor(n_bases/2))
    
    indexes_1 = np.arange(half)
    indexes_2 = np.arange(half, n_bases)

    return half, indexes_1, indexes_2


def make_pair_map(hb_id_1, hb_id_2):
    pair_map = {}
    for i, j in zip(hb_id_1, hb_id_2):
        pair_map[i] = j
        pair_map[j] = i
    return pair_map


def get_idx_of_coding_complement(indexes_1, indexes_2, pair_map, coding_seq_str):
    assert len(coding_seq_str) <= len(indexes_1), "coding sequence is longer than half of the structure"
    
    second_half_reversed = indexes_2[::-1]
    coding_indexes = second_half_reversed[:len(coding_seq_str)]
    coding_complement = [pair_map[i] for i in coding_indexes]
    
    return coding_indexes, coding_complement


def find_where_coding_binds_coding(coding_indexes, coding_complement, pair_map):
    coding_binds_coding_indexes = []
    for code in coding_indexes:
        if code in coding_complement:
            coding_binds_coding_indexes.append((code, pair_map[code]))
    return coding_binds_coding_indexes


def get_index_to_seq_map(strucutre):
    index_to_seq_map = {}
    for idx in range(strucutre.get_num_bases()):
        base = strucutre.get_base(idx)
        index_to_seq_map[idx] = base.base[0]
    return index_to_seq_map


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