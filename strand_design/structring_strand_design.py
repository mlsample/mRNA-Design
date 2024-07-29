import ipy_oxdna as iox
import ipy_oxdna.dna_structure as dna
from oxDNA_analysis_tools.output_bonds import output_bonds
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, inbox, describe
from oxDNA_analysis_tools.UTILS.oxview import oxdna_conf
from folder_base import get_energy_ratio, setup_ssorigami_from_files

import os
import multiprocessing as mp
import sys
import numpy as np
import argparse
import subprocess
from pathlib import Path
import pandas as pd
from copy import deepcopy
import json
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main():
    strucutre_file, topology_file, input_md_file, coding_sequence, output_name, traj_file, five_prime, three_prime = collect_args()

    coding_seq_str = parse_coding_sequence(coding_sequence)

    strucutre = get_structure(strucutre_file, topology_file)
    n_bases = len(strucutre.strands[0])

    hb_id_1, hb_id_2 = get_hblist(strucutre_file, topology_file, input_md_file, traj_file, n_bases)
    
    half, indexes_1, indexes_2 = get_halfs(strucutre)
    
    pair_map = make_pair_map(hb_id_1, hb_id_2)
    
    coding_indexes, coding_with_complement, coding_complement, coding_no_complement = get_idx_of_coding_complement(indexes_1, indexes_2, pair_map, coding_seq_str)
    
    coding_binds_coding_indexes = find_where_coding_binds_coding(coding_indexes, coding_complement, pair_map)
    
    mutated_struct = mutate_coding(strucutre, coding_indexes, coding_seq_str)
    
    index_to_seq_map_with_coding = get_index_to_seq_map(mutated_struct)
    
    mutate_noncoding(mutated_struct, index_to_seq_map_with_coding, coding_with_complement,
                     coding_complement, coding_no_complement, pair_map,
                     coding_binds_coding_indexes, coding_indexes)
    
    stats = validate_mutation(mutated_struct, index_to_seq_map_with_coding, coding_with_complement,
                      coding_complement, coding_no_complement, pair_map,
                      coding_binds_coding_indexes, coding_seq_str, coding_indexes, strucutre)
    
    ori_new_tocolor = plot_statistics(mutated_struct, pair_map, stats, mutated_struct, topology_file, output_name)
    
    add_5_prime(five_prime, mutated_struct)
    
    add_3_prime(three_prime, mutated_struct)
    
    export_structure(mutated_struct, strucutre, output_name, stats, ori_new_tocolor, five_prime, three_prime)
    
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
                        help='Path to .txt or .fasta file coding sequence with bases in 5` to 3` order EX: ./coding_seq.txt')
    parser.add_argument('-o', '--output_name', metavar='output_name', type=str,
                        help='name of the Default: ./fatcat_tmalign_homology_search.csv')
    parser.add_argument('--traj_file', metavar='traj_file', type=str,
                        help='Path to oxDNA .dat trajectory file of for determining hb pairs')
    parser.add_argument('--five_prime', metavar='five_prime', type=str,
                        help='Path to a .txt or .fasta file with bases in 5` to 3` order to add to 5 prime end')
    parser.add_argument('--three_prime', metavar='three_prime', type=str,
                        help='Path to a .txt or .fasta file with bases with bases in 5` to 3` order to add to 3 prime end')
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
        output_name = Path('coding_sequence_embbeded_strucutre')
    else:
        output_name = Path(args.output_name)
        
    
    if not args.traj_file:
        traj_file = None
    else:
        traj_file = Path(args.traj_file)
        if not traj_file.exists():
            raise ValueError('Trajectory file does not exist')
        
    if not args.five_prime:
        five_prime = None
    else:
        five_prime = Path(args.five_prime)
        if not five_prime.exists():
            raise ValueError('Five prime file does not exist')
    
    if not args.three_prime:
        three_prime = None
    else:
        three_prime = Path(args.three_prime)
        if not three_prime.exists():
            raise ValueError('Three prime file does not exist')
    
    if not args.force_overwrite:
        if (Path(f'{output_name}.dat').exists()) or (Path(f'{output_name}.top').exists()):
            raise ValueError(f'Output file {output_name} already exists.\
                Please use the --force_overwrite flag to overwrite the file.')
    return strucutre_file, topology_file, input_md_file, coding_sequence, output_name, traj_file, five_prime, three_prime


def parse_coding_sequence(coding_sequence):
    if coding_sequence.suffix == '.fasta':
        with open(coding_sequence, 'r') as f:
            lines = f.readlines()
            coding_seq_str = lines[1]
    elif coding_sequence.suffix == '.txt':
        with open(coding_sequence, 'r') as f:
            coding_seq_str = f.readline()

    unique_bases = set(coding_seq_str)
    possible_bases = set('ATCGUatcgu')

    assert unique_bases.issubset(possible_bases), "coding_sequence must only contain A, T, C, G, or U"
    coding_seq_str = coding_seq_str.upper()
    return coding_seq_str


def get_structure(strucutre_file, topology_file):
    strucutre = dna.load_dna_structure(topology_file, strucutre_file)
    return strucutre


def run_output_bonds(input_md_file, strucutre_file):
    p1 = f'oat output_bonds "{input_md_file.as_posix()}" "{strucutre_file.as_posix()}" '

    p2 = """| grep -v "#" | gawk '{if($7 < -0.1){print $1 " " $2 " " $7 " "}}' > hblist.txt"""

    invovation = p1 + p2
    start_dir = os.getcwd()
    os.chdir(input_md_file.parent)
    try:
        result = subprocess.run(invovation, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(e)
        raise ValueError('The run output bonds script failed')
    os.chdir(start_dir)
    if result.returncode != 0:
        print('The output_bonds script failed')
        sys.exit(1)
    
    
    return None

def get_hblist(strucutre_file, topology_file, input_md_file, traj_file, n_bases):
    if traj_file is not None:
        strucutre_file = traj_file

    # run_output_bonds(input_md_file, strucutre_file)
    with open('hblist.txt', 'r') as f:
        lines = f.readlines()
    lines_strip = [line.strip() for line in lines]
    lines_split = [line.split(' ') for line in lines_strip]
    lines_array = []
    for line in lines_split:
        try:
            lines_array.append(np.array(line, dtype=float))
        except:
            pass
            # print(line)
            # raise ValueError('The run output bonds script ')
    lines_one_array = np.array(lines_array)
    columns = ['id1', 'id2', 'HB']
    df = pd.DataFrame(lines_one_array, columns=columns, dtype=float)
    df_result = df.groupby(['id1', 'id2']).mean().reset_index()
    unq_counts_1 = np.unique(df_result['id1'], return_counts=True)
    unq_counts_2 = np.unique(df_result['id2'], return_counts=True)
    idx_1 = np.where(unq_counts_1[1] > 1)
    idx_2 = np.where(unq_counts_2[1] > 1)
    failed_1 = unq_counts_1[0][np.where(unq_counts_1[1] > 1)] 
    failed_2 = unq_counts_2[0][np.where(unq_counts_2[1] > 1)] 

    indexes_to_drop = []
    for fail in failed_1:
        rows_to_look_at = df_result[df_result['id1'] == fail]
        indexes_to_drop.append(rows_to_look_at[rows_to_look_at['HB'] != rows_to_look_at['HB'].min()].index)
    for fail in failed_2:
        rows_to_look_at = df_result[df_result['id2'] == fail]
        indexes_to_drop.append(rows_to_look_at[rows_to_look_at['HB'] != rows_to_look_at['HB'].min()].index)
         
    index_pd = np.unique(indexes_to_drop)
    df_result = df_result.drop(index=index_pd)
    df_result = df_result.reset_index(drop=True)
        
    hb_id_1 = df_result['id1'].reset_index(drop=True).map(int)
    hb_id_2 = df_result['id2'].reset_index(drop=True).map(int)
    
    df_result.to_csv('hb_list_traj.csv')
    
    with open('hb_list_traj.txt', 'w') as f:
        for id1, id2 in zip(hb_id_1, hb_id_2):
            f.write(f'{id1} {id2}\n')

    return hb_id_1, hb_id_2


def get_halfs(strucutre):
    n_bases = strucutre.get_num_bases()
    half = int(np.floor(n_bases/2))
    
    indexes_1 = np.array(np.arange(half), dtype=int)
    indexes_2 = np.array(np.arange(half, n_bases), dtype=int)

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
    
    coding_complement = []
    coding_indexes = second_half_reversed[:len(coding_seq_str)]
    
    coding_with_complement = [i for i in coding_indexes if int(i) in pair_map]
    coding_complement = [pair_map[int(i)] for i in coding_indexes if int(i) in pair_map]
    
    coding_no_complement = [i for i in coding_indexes if int(i) not in pair_map]
    
    return coding_indexes, coding_with_complement, coding_complement, coding_no_complement


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


def mutate_coding(strucutre, coding_indexes, coding_seq_str):
    mutated_struct = deepcopy(strucutre)
    mutated_struct.strands[0].mutate_sequence(coding_seq_str, min(coding_indexes), max(coding_indexes)+1) 
    coding_str_U_to_T =  ''.join([base if base != 'U' else 'T' for base in coding_seq_str])
    struct_bases = ''.join(mutated_struct.strands[0].bases[coding_indexes].tolist())
    assert coding_str_U_to_T == struct_bases
    return mutated_struct


def mutate_noncoding(mutated_struct, index_to_seq_map_with_coding, coding_with_complement, coding_complement, coding_no_complement, pair_map, coding_binds_coding_indexes, coding_indexes):
    # After I mutated the coding, I now want to mutate non-coding that is complement to the coding, but which is not itself coding
    i = 0
    j = 0
    coding_binds_coding = np.unique(list(map(int, [vals for tpl in coding_binds_coding_indexes for vals in tpl])))
    for code_w_comp, code_comp in zip(coding_with_complement, coding_complement):
        if (code_w_comp not in coding_binds_coding) and (code_comp not in coding_binds_coding):
            j += 1
            assert code_comp not in coding_indexes
            assert code_w_comp not in coding_no_complement
            assert pair_map[code_w_comp] == code_comp, 'Logic error'
            coding_base = index_to_seq_map_with_coding[code_w_comp]
            assert mutated_struct.strands[0].bases[code_w_comp] == coding_base
            new_base = get_compseq(coding_base)
            mutated_struct.strands[0].mutate_sequence(new_base, code_comp)
        else:
            i += 1
    
    return None


def validate_mutation(mutated_struct, index_to_seq_map_with_coding, coding_with_complement, coding_complement, coding_no_complement, pair_map, coding_binds_coding_indexes, coding_seq_str, coding_indexes, strucutre):
    mutated_struct_bak = deepcopy(mutated_struct)
    strucutre_bak = deepcopy(strucutre)
    
    struct_bases = mutated_struct_bak.strands[0].bases
    
    min_coding_idx = min(coding_indexes)
    max_coding_idx = max(coding_indexes)
    coding_slice = slice(min_coding_idx, max_coding_idx +1)
    
    struct_coding_bases = struct_bases[coding_slice]
    
    struct_obj_coding_bases = mutated_struct_bak.strands[0][coding_slice]

    real_coding_5_to_3 = ''.join([base if base != 'U' else 'T' for base in coding_seq_str])
    real_coding_3_to_5 = ''.join([seq if seq != 'U' else 'T' for seq in coding_seq_str][::-1])
    struct_bases_str_5_3 = ''.join(struct_coding_bases.tolist()[::-1])
    struct_bases_str_3_5 = ''.join(struct_coding_bases.tolist())
    
    assert real_coding_5_to_3 == struct_bases_str_5_3
    assert real_coding_3_to_5 == struct_bases_str_3_5
    
    coding_with_complement_bases = [mutated_struct_bak.strands[0][int(code_w_comp)] for code_w_comp in coding_with_complement]
    coding_w_comp_complements = [pair_map[int(code_w_comp)] for code_w_comp in coding_with_complement]
    coding_w_comp_complements_bases = [mutated_struct_bak.strands[0][int(comps)] for comps in coding_w_comp_complements]
    
    c_w_c_b_str = ''.join([base_obj.base for base_obj in coding_with_complement_bases])
    c_w_c_c_b_str = ''.join([base_obj.base for base_obj in coding_w_comp_complements_bases])
    
    c_w_c_c_b_str_comp = ''.join([get_compseq(base) for base in c_w_c_c_b_str])
    
    assert c_w_c_b_str != c_w_c_c_b_str_comp
    
    all_bad_idxes = [idx for tpl in coding_binds_coding_indexes for idx in tpl]
    
    coding_with_comp_not_self = [code_w_comp for code_w_comp in coding_with_complement if code_w_comp not in all_bad_idxes]
    the_comp_idxes = [pair_map[int(c_w_c_n_s)] for c_w_c_n_s in coding_with_comp_not_self]
    
    coding_with_comp_not_self_bases = [mutated_struct_bak.strands[0][int(base_idx)] for base_idx in coding_with_comp_not_self]
    the_comp_idxes_bases = [mutated_struct_bak.strands[0][int(base_idx)] for base_idx in the_comp_idxes]
    
    coding_with_comp_not_self_bases_str = ''.join([base_obj.base for base_obj in coding_with_comp_not_self_bases])
    the_comp_idxes_bases_str = ''.join([base_obj.base for base_obj in the_comp_idxes_bases])
    the_comp_idxes_bases_str_comp = ''.join([get_compseq(base) for base in the_comp_idxes_bases_str])
    
    assert coding_with_comp_not_self_bases_str == the_comp_idxes_bases_str_comp
    
    for key,value in pair_map.items():
        if (key not in all_bad_idxes) and (value not in all_bad_idxes):
            key_base = strucutre_bak.strands[0][int(key)].base
            value_base = strucutre_bak.strands[0][int(value)].base
            comp_value_base = get_compseq(value_base)
            assert comp_value_base == key_base, f'{key}, {value} failed'
        
    base_idx_set = set(range(len(mutated_struct_bak.strands[0])))
    paired_idx_set = set(pair_map.keys())
    unpaired_nucs = base_idx_set.difference(paired_idx_set)
    unpaired_in_coding = set(coding_indexes).intersection(unpaired_nucs)
    unpaired_in_coding_list = sorted(list(unpaired_in_coding))
    unpaired_in_strucutring = unpaired_nucs.difference(unpaired_in_coding)
    unpaired_in_strucutring_list = sorted(list(unpaired_in_strucutring))
    
    all_coding = set(map(int, coding_indexes))
    coding_no_comp = set(map(int, coding_no_complement))
    coding_self_binding = set(map(int, all_bad_idxes))
    coding_bound = all_coding.difference(coding_no_comp).difference(coding_self_binding)
    
    all_structuring = base_idx_set.difference(all_coding)
    structuring_no_comp = all_structuring.difference(paired_idx_set)
    structuring_bound_to_coding = set([pair_map[int(bound_coding)] for bound_coding in list(coding_bound)])
    structuring_selfbound = all_structuring.difference(structuring_bound_to_coding).difference(structuring_no_comp)
    
    assert len(set(coding_indexes).intersection(unpaired_nucs).difference(set(coding_no_complement))) == 0
    
    return coding_no_comp, coding_self_binding, coding_bound, structuring_no_comp, structuring_bound_to_coding, structuring_selfbound


def plot_statistics(mutated_strcutre, pair_map, stats, mutated_struct, topology_file, output_name):
    coding_no_comp, coding_self_binding, coding_bound, structuring_no_comp, structuring_bound_to_coding, structuring_selfbound = stats
    
    total_nucs = len(mutated_struct.strands[0])
    total_paired_nucs = len(pair_map)
    total_pairs = len(pair_map) // 2
    
    coding_total_unpaired = len(coding_self_binding)
    total_missing_pairs = coding_total_unpaired // 2
    
    total_no_comp = len(coding_no_comp) + len(structuring_no_comp)
    
    lost_pairs = total_no_comp - coding_total_unpaired
    

    hb_list = 'hb_list_traj.txt'
    
    
    top_file_path = output_name.with_suffix('.top') 
    dat_file_path = output_name.with_suffix('.dat')
    mutated_strcutre.export_top_conf(top_file_path, dat_file_path)

    original = setup_ssorigami_from_files(topology_file, hb_list)
    new = setup_ssorigami_from_files(output_name.with_suffix('.top'), hb_list)
    
    ori_new_diffs, ori_new_tocolor = get_energy_ratio(original, new, color_threshold=-0.4)

    label = f'Number of domains:{len(ori_new_diffs)}\nMean percent energy change: {np.mean(ori_new_diffs):.2f}'
    sns.histplot(ori_new_diffs, label=label)
    plt.xlabel('Percent energy change between original and new duplex domains')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_name.with_suffix('.png'))
    
    stat_print = f"""
Single-stranded origami strucutring strand design algorithm result info:
    Base struct nucs: {total_nucs}
    Base struct HB paired nucs: {total_paired_nucs}
    Base struct unpaired nucs: {total_nucs - total_paired_nucs}
    Base struct % HB paired = {100*(total_paired_nucs / total_nucs):.2f}%
    
    Coding struct HB paired nucs: {total_paired_nucs - coding_total_unpaired} 
    Coding struct HB lost paired nucs: {coding_total_unpaired}
    Coding struct % HB paired: {100*((total_paired_nucs - coding_total_unpaired) / total_nucs):.2f}%
    
    Calculated domains: {len(ori_new_diffs)}
    Destabilized domains: {len([ori for ori in ori_new_diffs if ori < -0.4])}
    Mean percent domain energy change including unchanged domains: {np.mean(ori_new_diffs)*100:.2f} %
    Mean percent domain energy change excluding unchanged domains: {np.mean([ori for ori in ori_new_diffs if ori != 0])*100:.2f} %
    
    Percent domain energy change histogram saved to: {output_name.with_suffix('.png')}
    """
    
    print(stat_print)
    print(f'    Info saved to: {output_name.with_suffix(".txt")}\n')
    
    with open(output_name.with_suffix('.txt'), 'w') as f:
        f.write(stat_print)
    
    return ori_new_tocolor


def add_5_prime(five_prime, mutated_strcutre):
    if five_prime is not None:
        current_five_prime_base = mutated_strcutre.strands[0][-1]

        bases_to_add = parse_coding_sequence(five_prime)
        bases_to_add = ''.join([base if base != 'U' else 'T' for base in bases_to_add][::-1])
        fwd_strand, reverse_strand = dna.construct_strands(bases_to_add, current_five_prime_base.pos - 0.3897, current_five_prime_base.a3)

        mutated_strcutre.strands[0].append(fwd_strand)
    return None


def add_3_prime(three_prime, mutated_strcutre):
    if three_prime is not None:
        current_three_prime_base = mutated_strcutre.strands[0][0]

        bases_to_add = parse_coding_sequence(three_prime)
        bases_to_add = ''.join([base if base != 'U' else 'T' for base in bases_to_add][::-1])
        coord = current_three_prime_base.a3
        pos_normed = current_three_prime_base.pos / np.linalg.norm(current_three_prime_base.pos)
        contribution =  current_three_prime_base.pos - coord * (0.3897 * len(bases_to_add))
        
        fwd_strand, reverse_strand = dna.construct_strands(bases_to_add, contribution, coord)

        mutated_strcutre.strands[0].prepend(fwd_strand)
    return None


def export_structure(mutated_strcutre, strucutre, output_name, stats, ori_new_tocolor, five_prime, three_prime):
    top_file_path = output_name.with_suffix('.top') 
    dat_file_path = output_name.with_suffix('.dat')
    mutated_strcutre.export_top_conf(top_file_path, dat_file_path)
    
    mutated_strcutre.export_oxview(output_name.with_suffix('.oxview'))
    with open(output_name.with_suffix('.oxview'), 'r') as f:
        oxview_file = json.load(f)
    
    coding_no_comp, coding_self_binding, coding_bound, structuring_no_comp, structuring_bound_to_coding, structuring_selfbound = stats
    
    monomers  = oxview_file['systems'][0]['strands'][0]['monomers']
    selections = oxview_file['selections']
    
    # Define a list of color names
    color_names = [
        "lightcoral", "darkturquoise", "mediumpurple", "maroon", "orchid", "plum", "snow", "steelblue", "lightsteelblue"
    ]
    
    # Convert color names to RGB integers
    color_pallete = [color_name_to_rgb(color) for color in color_names]
    
    regions_to_color = [coding_self_binding, coding_bound, structuring_bound_to_coding, ori_new_tocolor, structuring_no_comp, coding_no_comp, structuring_selfbound]
    regions_to_color = [np.array(list(indexes)) for indexes in regions_to_color]
    
    if three_prime is not None:
        three_prime_bases_added = len(parse_coding_sequence(three_prime))
        for indexes in regions_to_color:
            indexes += three_prime_bases_added
        regions_to_color.append(list(range(three_prime_bases_added)))
    
    if five_prime is not None:
        bases_added = len(parse_coding_sequence(five_prime))
        if three_prime is not None:
            starting_index = three_prime_bases_added + len(strucutre.strands[0])
        else:
            starting_index = len(strucutre.strands[0])
        five_idxes = list(range(starting_index, starting_index + bases_added))    
        regions_to_color.append(five_idxes)
            
    for color_idx, indexes in enumerate(regions_to_color):
        my_color = color_pallete[color_idx]
        for idx in indexes:
            monomers[idx]['color'] = my_color 
        
    sys_names = ["coding_self_binding", "coding_bound_to_structuring", "structuring_bound_to_coding", "highly_destabilized", "structuring_no_comp", "coding_no_comp", "structuring_self_bound"]
    if three_prime is not None:
        sys_names.append('three_prime')
    if five_prime is not None:
        sys_names.append('five_prime')
    
    for sys_idx, indexes in enumerate(regions_to_color):
        select = sys_names[sys_idx].replace('_', ' ')
        select += f": {color_names[sys_idx]}"
        my_select = [select, list(map(int, list(indexes)))]
        selections.append(my_select)
        
    oxview_file['systems'][0]['strands'][0]['monomers'] = monomers
    oxview_file['selections'] = selections
    
    with open(output_name.with_suffix('.oxview'), 'w') as f:
        json.dump(oxview_file, f)
    
    
    struct_bases = mutated_strcutre.strands[0].bases
    struct_bases_str_5_3 = ''.join(struct_bases.tolist()[::-1])
    
    regions_to_color_mirrored = mirror_indexes(regions_to_color, len(struct_bases))
    annotations = []                                    
    for color_idx, indexes in enumerate(regions_to_color_mirrored):
        my_name = sys_names[color_idx]
        grouped = group_consecutive_indexes(indexes)
        for sub_idxes in grouped:
        
            min_idx = int(min(sub_idxes))
            max_idx = int(max(sub_idxes)) +1
            annot = (min_idx, max_idx, my_name, 'DNA')
            annotations.append(annot)
            
    coding = np.concatenate([regions_to_color_mirrored[0], regions_to_color_mirrored[1], regions_to_color_mirrored[5]])
    structring = np.concatenate([regions_to_color_mirrored[2], regions_to_color_mirrored[4], regions_to_color_mirrored[6]])
        
    annotations.append((int(min(coding)), int(max(coding))+1, 'all_coding', 'DNA'))
    annotations.append((int(min(structring)), int(max(structring))+1, 'all_structuring', 'DNA'))
    
    create_benchling_dna_file(struct_bases_str_5_3, annotations, output_name.with_suffix('.dna'), output_name)
    
    return None


def create_benchling_dna_file(sequence, annotations, file_name, output_name):
    """
    Create a Benchling-readable .dna file.
    
    :param sequence: str, the full DNA sequence.
    :param annotations: list of tuples, each containing (start, end, label, type).
    :param file_name: str, the name of the output .dna file.
    """
    # Create a SeqRecord object
    seq_record = SeqRecord(Seq(sequence), id="ExampleID", name=f"{output_name.stem}", description=f"Annotated {output_name}")
    seq_record.annotations["molecule_type"] = "DNA"
    # Add annotations
    for start, end, label, feature_type in annotations:
        feature = SeqFeature(FeatureLocation(start, end), type=feature_type, qualifiers={"label": label})
        seq_record.features.append(feature)
    
    # Write to file
    SeqIO.write(seq_record, file_name, "genbank")


def mirror_indexes(index_lists, max_index):
    mirrored_lists = []
    for lst in index_lists:
        mirrored_list = [max_index - index for index in lst]
        mirrored_lists.append(mirrored_list)
    return mirrored_lists


def group_consecutive_indexes(indexes):
    if not indexes:
        return []

    # Sort the indexes if they are not sorted
    indexes.sort()

    # Initialize variables
    grouped_indexes = []
    current_group = [indexes[0]]

    # Iterate through the list starting from the second element
    for i in range(1, len(indexes)):
        if indexes[i] == indexes[i-1] + 1:
            # Current index is consecutive, add to the current group
            current_group.append(indexes[i])
        else:
            # Current index is not consecutive, save the current group and start a new one
            grouped_indexes.append(current_group)
            current_group = [indexes[i]]

    # Add the last group to the result
    grouped_indexes.append(current_group)

    return grouped_indexes

def rgb_to_int(r, g, b):
    return (r << 16) + (g << 8) + b


def color_name_to_rgb(color_name):
    rgb = mcolors.to_rgb(color_name)
    return rgb_to_int(*(int(x * 255) for x in rgb))


def get_compseq(seq):
    complement_translation = str.maketrans('ACGT', 'TGCA')
    return seq.translate(complement_translation)[::-1]


if __name__ == '__main__':
    main()