#!/usr/bin/env python3
"""
Slice the structuring half of a coding ssOrigami into DNA origami staples.

The input is the structure written by ``strand_design/structring_strand_design.py``: a single
strand whose 5` half carries the coding sequence (optionally flanked by UTRs) and whose 3` half
is the structuring complement that folds the nanostructure.

oxDNA stores the nucleotides of a strand 3` -> 5`, so index 0 is the 3` terminus.  In that index
space the design script's output is laid out as::

    index   0        L3-1  L3                 C0-1  C0           C1-1  C1        N-1
            |  3` UTR  |  |  structuring region  |  |  coding sequence  |  |  5` UTR  |
    5` -> 3` reading order is the reverse: 5` UTR, coding, structuring, 3` UTR

This script nicks that strand into:

* one **coding strand** -- the 5` UTR + coding sequence, with the 3` UTR moved from the far end
  of the structuring region to the coding sequence`s own 3` end (where the coding/structuring
  nick is made) and re-oriented to continue off that terminus.
* many **staples** -- the structuring region cut every ``--staple_length`` nucleotides.

Several slicing methods are supported (``--method``); ``uniform`` is the simple fixed-length cut.

Outputs (for ``--output_name NAME``):

    NAME.top, NAME.dat    the sliced structure in oxDNA format
    NAME.oxview           oxView file coloured by strand with a selections legend
    NAME_strands.csv      primary sequences (5` -> 3`) of every strand, with length/GC/Tm
    NAME_strands.fasta    the same sequences as FASTA
    NAME_coding.dna       Benchling-readable (GenBank) coding strand with UTR annotations
    NAME.txt              summary statistics
"""

import argparse
import csv
import json
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import matplotlib.colors as mcolors
import ipy_oxdna.strucutre_editor.dna_structure as dna
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt

SLICING_METHODS = ('uniform', 'even')

# Colours used for the oxView export, one per region / cycled over the staples
REGION_COLORS = {'five_prime': 'lightsteelblue', 'coding': 'mediumpurple', 'three_prime': 'plum'}
STAPLE_COLORS = ('lightcoral', 'darkturquoise', 'steelblue', 'orchid')


@dataclass
class Fragment:
    """One strand of the sliced structure, with a record of where it came from."""
    name: str
    kind: str                                   # 'coding', 'staple' or 'three_prime'
    strand: dna.DNAStructureStrand
    source_start: int                           # inclusive, index in the input single strand
    source_stop: int                            # exclusive
    note: str = ''
    sub_regions: list = field(default_factory=list)   # (start, stop, label) in local strand index


def main():
    args = collect_args()

    seqs = parse_region_sequences(args)

    strucutre = dna.load_dna_structure(args.topology_file, args.strucutre_file)
    assert strucutre.get_num_strands() == 1, \
        f'Expected a single stranded origami, got {strucutre.get_num_strands()} strands'
    strand = strucutre.strands[0]

    regions = get_regions(len(strand), seqs)
    verify_regions(strand, regions, seqs, args)

    fragments = slice_strand(strand, regions, args)

    validate_slicing(strand, regions, fragments, args)

    sliced = dna.DNAStructure([frag.strand for frag in fragments],
                              strucutre.time, strucutre.box, strucutre.energy)

    export_structure(sliced, fragments, args)
    export_sequences(fragments, args)
    export_coding_annotation(fragments, args)
    write_summary(strand, regions, fragments, args)

    return 0


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Nicks the structuring half of a coding ssOrigami into DNA origami staples',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--strucutre_file', metavar='strucutre_file', type=str,
                        help='Path to the oxDNA .dat file written by structring_strand_design.py')
    parser.add_argument('-t', '--topology_file', metavar='topology_file', type=str,
                        help='Path to the oxDNA .top file written by structring_strand_design.py')
    parser.add_argument('-c', '--coding_sequence', metavar='coding_sequence', type=str,
                        help='Path to the .txt or .fasta coding sequence used by the design script')
    parser.add_argument('--five_prime', metavar='five_prime', type=str,
                        help='Path to the .txt or .fasta 5 prime sequence used by the design script')
    parser.add_argument('--three_prime', metavar='three_prime', type=str,
                        help='Path to the .txt or .fasta 3 prime sequence used by the design script')
    parser.add_argument('--coding_length', metavar='coding_length', type=int,
                        help='Length of the coding region, instead of passing --coding_sequence')
    parser.add_argument('--five_prime_length', metavar='five_prime_length', type=int,
                        help='Length of the 5 prime region, instead of passing --five_prime')
    parser.add_argument('--three_prime_length', metavar='three_prime_length', type=int,
                        help='Length of the 3 prime region, instead of passing --three_prime')
    parser.add_argument('-n', '--staple_length', metavar='staple_length', type=int, default=35,
                        help='Number of nucleotides between nicks in the structuring region')
    parser.add_argument('--min_staple_length', metavar='min_staple_length', type=int, default=15,
                        help='Shortest staple allowed; a shorter leftover is merged into its neighbour')
    parser.add_argument('-m', '--method', metavar='method', type=str, default='uniform',
                        choices=SLICING_METHODS,
                        help='How to place the nicks: "uniform" cuts exactly every --staple_length '
                             'nucleotides, "even" spreads the region over staples of near-equal length')
    parser.add_argument('--nick_from', metavar='nick_from', type=str, default='coding_junction',
                        choices=('coding_junction', 'three_prime_junction'),
                        help='Which end of the structuring region the fixed-length cuts start from')
    parser.add_argument('--three_prime_placement', metavar='three_prime_placement', type=str,
                        default='coding_3p', choices=('coding_3p', 'keep'),
                        help='"coding_3p" moves the 3 prime region onto the 3 prime end of the coding '
                             'strand, "keep" leaves it as its own strand where the design script put it')
    parser.add_argument('--coding_as_rna', action='store_true',
                        help='Report the coding strand sequence as RNA (T -> U) in the sequence outputs')
    parser.add_argument('--no_verify_regions', action='store_true',
                        help='Skip checking that the structure bases match the provided sequence files')
    parser.add_argument('-o', '--output_name', metavar='output_name', type=str,
                        default='structuring_staples',
                        help='Basename of the output files')
    parser.add_argument('--force_overwrite', action='store_true',
                        help='Force overwrite of existing results if they exist.')

    args = parser.parse_args()
    return args, parser


def collect_args():
    """
    Checks that the provided arguments exist and are self consistent.

    Raises:
        ValueError: if a required argument is missing or a provided path does not exist

    Returns:
        The parsed arguments, with every file argument turned into a Path
    """
    args, parser = parse_arguments()

    if not args.strucutre_file:
        parser.print_help()
        raise ValueError('Please provide a strucutre_file using the -s or --strucutre_file flag')
    args.strucutre_file = Path(args.strucutre_file)
    if not args.strucutre_file.exists():
        raise ValueError(f'Provided strucutre_file path {args.strucutre_file} does not exist')
    assert args.strucutre_file.suffix == '.dat', \
        f'strucutre_file must be a .dat file and it is currently {args.strucutre_file.suffix}'

    if not args.topology_file:
        parser.print_help()
        raise ValueError('Please provide a topology_file using the -t or --topology_file flag')
    args.topology_file = Path(args.topology_file)
    if not args.topology_file.exists():
        raise ValueError(f'Provided topology_file path {args.topology_file} does not exist')
    assert args.topology_file.suffix == '.top', 'topology_file must be a .top file'

    for flag, length_flag in (('coding_sequence', 'coding_length'),
                              ('five_prime', 'five_prime_length'),
                              ('three_prime', 'three_prime_length')):
        path = getattr(args, flag)
        if path is None:
            continue
        path = Path(path)
        if not path.exists():
            raise ValueError(f'Provided {flag} path {path} does not exist')
        if path.suffix not in ['.txt', '.fasta']:
            raise ValueError(f'{flag} must be a .txt or .fasta file')
        if getattr(args, length_flag) is not None:
            raise ValueError(f'Provide either --{flag} or --{length_flag}, not both')
        setattr(args, flag, path)

    if args.coding_sequence is None and args.coding_length is None:
        parser.print_help()
        raise ValueError('Please provide the coding region using --coding_sequence or --coding_length')

    if args.staple_length < 1:
        raise ValueError('--staple_length must be at least 1')
    if args.min_staple_length < 1:
        raise ValueError('--min_staple_length must be at least 1')
    if args.min_staple_length > args.staple_length:
        raise ValueError('--min_staple_length cannot be larger than --staple_length')

    args.output_name = Path(args.output_name)
    if not args.force_overwrite:
        existing = [path for path in output_paths(args.output_name) if path.exists()]
        if existing:
            raise ValueError(f'Output file(s) {", ".join(map(str, existing))} already exist. '
                             'Please use the --force_overwrite flag to overwrite them.')

    return args


def output_paths(output_name):
    """Every file this script writes, for the overwrite check and the summary."""
    return [output_name.with_suffix('.top'),
            output_name.with_suffix('.dat'),
            output_name.with_suffix('.oxview'),
            output_name.with_suffix('.txt'),
            output_name.parent / f'{output_name.stem}_strands.csv',
            output_name.parent / f'{output_name.stem}_strands.fasta',
            output_name.parent / f'{output_name.stem}_coding.dna']


def parse_sequence_file(sequence_file):
    """
    Reads a one line .txt or a two line .fasta of bases in 5` to 3` order.

    This mirrors parse_coding_sequence in strand_design/structring_strand_design.py, but it also
    strips trailing whitespace so the length can be compared against the structure.
    """
    lines = [line.strip() for line in sequence_file.read_text().splitlines() if line.strip()]
    if sequence_file.suffix == '.fasta':
        assert len(lines) >= 2, 'fasta file must have a header line and a sequence line'
        sequence = ''.join(lines[1:])
    else:
        assert len(lines) == 1, 'txt file must have 1 line'
        sequence = lines[0]

    assert set(sequence).issubset(set('ATCGUatcgu')), \
        f'{sequence_file} must only contain A, T, C, G, or U'
    return sequence.upper()


def to_dna(sequence):
    return sequence.upper().replace('U', 'T')


def to_rna(sequence):
    return sequence.upper().replace('T', 'U')


def parse_region_sequences(args):
    """Returns the coding / 5` / 3` sequences (as DNA), using '' for regions that are absent."""
    seqs = {}
    for name, length_flag in (('coding', 'coding_length'),
                              ('five_prime', 'five_prime_length'),
                              ('three_prime', 'three_prime_length')):
        path = getattr(args, 'coding_sequence' if name == 'coding' else name)
        length = getattr(args, length_flag)
        if path is not None:
            seqs[name] = to_dna(parse_sequence_file(path))
        elif length is not None:
            if length < 0:
                raise ValueError(f'--{length_flag} cannot be negative')
            seqs[name] = 'N' * length          # length only, contents unknown
        else:
            seqs[name] = ''
    return seqs


def get_regions(n_bases, seqs):
    """
    Maps the design script's layout onto index ranges of the single strand.

    Every range is a half open (start, stop) in the strand's own 3` -> 5` index space.
    """
    len_three = len(seqs['three_prime'])
    len_five = len(seqs['five_prime'])
    len_coding = len(seqs['coding'])

    coding_start = n_bases - len_five - len_coding
    if coding_start < len_three:
        raise ValueError(f'The 3 prime ({len_three}), coding ({len_coding}) and 5 prime ({len_five}) '
                         f'regions do not fit in the {n_bases} nucleotide structure')

    regions = {'three_prime': (0, len_three),
               'structuring': (len_three, coding_start),
               'coding': (coding_start, coding_start + len_coding),
               'five_prime': (n_bases - len_five, n_bases)}

    if regions['structuring'][1] <= regions['structuring'][0]:
        raise ValueError('The structuring region is empty, there is nothing to slice')

    return regions


def region_seq_5to3(strand, region):
    """The sequence of a region read 5` -> 3`, i.e. against the stored index order."""
    start, stop = region
    return ''.join(strand.bases[start:stop].tolist()[::-1])


def verify_regions(strand, regions, seqs, args):
    """Checks that the region boundaries really do line up with the sequences that were designed in."""
    if args.no_verify_regions:
        return None

    for name in ('five_prime', 'coding', 'three_prime'):
        expected = seqs[name]
        if not expected or set(expected) == {'N'}:
            continue                                    # length was given, not a sequence
        found = region_seq_5to3(strand, regions[name])
        if found != expected:
            mismatches = [i for i, (a, b) in enumerate(zip(expected, found)) if a != b]
            raise ValueError(
                f'The {name} region of {args.topology_file} does not match the {name} sequence '
                f'that was provided ({len(mismatches)} mismatching bases, first at position '
                f'{mismatches[0] if mismatches else "n/a"} of {len(expected)}). '
                'Check that the sequence files are the ones used to build this structure, or pass '
                '--no_verify_regions to skip this check.')
    return None


def chunk_lengths(total, staple_length, method, min_staple_length):
    """
    Splits ``total`` nucleotides into staple lengths.

    'uniform' cuts every ``staple_length`` nucleotides and merges a too short leftover into its
    neighbour, 'even' spreads the region over staples that differ by at most one nucleotide.
    """
    if total <= staple_length:
        return [total]

    if method == 'uniform':
        lengths = [staple_length] * (total // staple_length)
        remainder = total - sum(lengths)
        if remainder:
            if remainder < min_staple_length:
                lengths[-1] += remainder
            else:
                lengths.append(remainder)
    elif method == 'even':
        n_staples = max(1, round(total / staple_length))
        base, remainder = divmod(total, n_staples)
        lengths = [base + 1] * remainder + [base] * (n_staples - remainder)
    else:
        raise ValueError(f'Unknown slicing method {method}')

    assert sum(lengths) == total, 'Logic error, the staple lengths do not add up'
    return lengths


def structuring_bounds(regions, args):
    """
    Returns the (start, stop) index ranges of the staples, ordered 5` -> 3` along the
    structuring region, so that the first staple is the one at the coding junction.
    """
    start, stop = regions['structuring']
    lengths = chunk_lengths(stop - start, args.staple_length, args.method, args.min_staple_length)

    bounds = []
    if args.nick_from == 'coding_junction':
        cursor = stop
        for length in lengths:
            bounds.append((cursor - length, cursor))
            cursor -= length
        assert cursor == start, 'Logic error, the staples do not cover the structuring region'
    else:
        cursor = start
        for length in lengths:
            bounds.append((cursor, cursor + length))
            cursor += length
        assert cursor == stop, 'Logic error, the staples do not cover the structuring region'
        bounds = bounds[::-1]

    return bounds


def graft_three_prime(coding_strand, bases_3to5):
    """
    Rebuilds the 3` region as a straight helix continuing off the coding strand's 3` terminus.

    ``bases_3to5`` is the 3` region in the strand's own index order.  The helix is phased so the
    nucleotide that ends up adjacent to the junction continues the junction nucleotide's twist,
    and it is placed on the helical axis of the junction nucleotide rather than on its centre of
    mass (construct_strands offsets each nucleotide from the axis by CM_CENTER_DS along a1).
    """
    junction = coding_strand[0]

    helix_dir = np.array(junction.a3, dtype=float)
    helix_dir /= np.linalg.norm(helix_dir)

    a1 = np.array(junction.a1, dtype=float)
    axis_point = np.array(junction.pos, dtype=float) + dna.CM_CENTER_DS * a1
    start_pos = axis_point - helix_dir * (dna.BASE_BASE * len(bases_3to5))

    perp = a1 - helix_dir * np.dot(helix_dir, a1)
    perp /= np.linalg.norm(perp)
    # construct_strands turns a1 by one base pair per nucleotide going 5`, so wind the starting
    # a1 back by one turn per added nucleotide to land on the junction's phase
    step = dna.get_rotation_matrix(helix_dir.copy(), 1, 'bp')
    perp = np.linalg.matrix_power(step.T, len(bases_3to5)) @ perp

    fwd_strand, _ = dna.construct_strands(bases_3to5, start_pos, helix_dir.copy(), perp=perp)
    coding_strand.prepend(fwd_strand)
    return None


def slice_strand(strand, regions, args):
    """Nicks the single strand into the coding strand and the structuring staples."""
    fragments = []

    # The coding sequence and the 5` region are contiguous, so they survive the nicking as one strand
    coding_start = regions['coding'][0]
    coding_stop = regions['five_prime'][1]
    coding_strand = deepcopy(strand[coding_start:coding_stop])

    len_coding = regions['coding'][1] - regions['coding'][0]
    len_five = regions['five_prime'][1] - regions['five_prime'][0]
    len_three = regions['three_prime'][1] - regions['three_prime'][0]

    moved_three_prime = len_three > 0 and args.three_prime_placement == 'coding_3p'
    if moved_three_prime:
        three_prime_bases = ''.join(strand.bases[slice(*regions['three_prime'])].tolist())
        graft_three_prime(coding_strand, three_prime_bases)

    offset = len_three if moved_three_prime else 0
    sub_regions = [(offset, offset + len_coding, 'coding'),
                   (offset + len_coding, offset + len_coding + len_five, 'five_prime')]
    if moved_three_prime:
        sub_regions.insert(0, (0, len_three, 'three_prime'))

    fragments.append(Fragment(
        name=f'{args.output_name.stem}_coding',
        kind='coding',
        strand=coding_strand,
        source_start=coding_start,
        source_stop=coding_stop,
        note='3 prime region moved to the coding 3 prime end' if moved_three_prime else '',
        sub_regions=sub_regions))

    bounds = structuring_bounds(regions, args)
    width = len(str(len(bounds)))
    for number, (start, stop) in enumerate(bounds, 1):
        fragments.append(Fragment(
            name=f'{args.output_name.stem}_staple_{number:0{width}d}',
            kind='staple',
            strand=deepcopy(strand[start:stop]),
            source_start=start,
            source_stop=stop))

    if len_three > 0 and not moved_three_prime:
        fragments.append(Fragment(
            name=f'{args.output_name.stem}_three_prime',
            kind='three_prime',
            strand=deepcopy(strand[slice(*regions['three_prime'])]),
            source_start=regions['three_prime'][0],
            source_stop=regions['three_prime'][1],
            note='3 prime region left where the design script put it'))

    return fragments


def validate_slicing(strand, regions, fragments, args):
    """Checks that nicking conserved every nucleotide and every sequence."""
    total_sliced = sum(len(frag.strand) for frag in fragments)
    assert total_sliced == len(strand), \
        f'Sliced structure has {total_sliced} nucleotides but the input had {len(strand)}'

    coding_fragment = fragments[0]
    expected_coding = (region_seq_5to3(strand, regions['five_prime'])
                       + region_seq_5to3(strand, regions['coding']))
    if args.three_prime_placement == 'coding_3p':
        expected_coding += region_seq_5to3(strand, regions['three_prime'])
    assert coding_fragment.strand.seq(from_5p=True) == expected_coding, \
        'The coding strand sequence changed while slicing'

    staples = [frag for frag in fragments if frag.kind == 'staple'
               and frag.source_start >= regions['structuring'][0]
               and frag.source_stop <= regions['structuring'][1]]
    stapled = ''.join(frag.strand.seq(from_5p=True) for frag in staples)
    assert stapled == region_seq_5to3(strand, regions['structuring']), \
        'The staples do not spell out the structuring region'

    covered = sorted((frag.source_start, frag.source_stop) for frag in staples)
    for (_, previous_stop), (start, _) in zip(covered, covered[1:]):
        assert previous_stop == start, 'The staples overlap or leave a gap'

    return None


def rgb_to_int(r, g, b):
    return (r << 16) + (g << 8) + b


def color_name_to_rgb(color_name):
    rgb = mcolors.to_rgb(color_name)
    return rgb_to_int(*(int(x * 255) for x in rgb))


def export_structure(sliced, fragments, args):
    """Writes the .top/.dat pair and an oxView file coloured by region and by staple."""
    top_file_path = args.output_name.with_suffix('.top')
    dat_file_path = args.output_name.with_suffix('.dat')
    sliced.export_top_conf(top_file_path, dat_file_path)

    oxview_path = args.output_name.with_suffix('.oxview')
    sliced.export_oxview(oxview_path)
    with open(oxview_path, 'r') as f:
        oxview_file = json.load(f)

    strands_json = oxview_file['systems'][0]['strands']
    assert len(strands_json) == len(fragments), 'oxView export does not match the sliced fragments'

    selections = oxview_file['selections']
    staple_number = 0
    for frag, strand_json in zip(fragments, strands_json):
        monomers = strand_json['monomers']
        if frag.kind == 'coding':
            for start, stop, label in frag.sub_regions:
                color_name = REGION_COLORS[label]
                ids = [monomers[idx]['id'] for idx in range(start, stop)]
                for idx in range(start, stop):
                    monomers[idx]['color'] = color_name_to_rgb(color_name)
                selections.append([f'{label.replace("_", " ")}: {color_name}', ids])
        else:
            if frag.kind == 'three_prime':
                color_name = REGION_COLORS['three_prime']
            else:
                color_name = STAPLE_COLORS[staple_number % len(STAPLE_COLORS)]
                staple_number += 1
            for monomer in monomers:
                monomer['color'] = color_name_to_rgb(color_name)
            selections.append([f'{frag.name}: {color_name}',
                               [monomer['id'] for monomer in monomers]])

    oxview_file['selections'] = selections
    with open(oxview_path, 'w') as f:
        json.dump(oxview_file, f)

    return None


def fragment_sequence(frag, args):
    """The primary sequence of a fragment, 5` -> 3`."""
    sequence = frag.strand.seq(from_5p=True)
    if frag.kind == 'coding' and args.coding_as_rna:
        sequence = to_rna(sequence)
    return sequence


def gc_percent(sequence):
    bases = sequence.upper()
    if not bases:
        return 0.0
    return 100.0 * (bases.count('G') + bases.count('C')) / len(bases)


def melting_temperature(sequence):
    try:
        return float(mt.Tm_NN(Seq(to_dna(sequence))))
    except (ValueError, KeyError, IndexError):
        return float('nan')


def export_sequences(fragments, args):
    """Writes the primary sequences of every strand as a CSV and as a FASTA."""
    csv_path = args.output_name.parent / f'{args.output_name.stem}_strands.csv'
    fasta_path = args.output_name.parent / f'{args.output_name.stem}_strands.fasta'

    rows = []
    for frag in fragments:
        sequence = fragment_sequence(frag, args)
        rows.append({'name': frag.name,
                     'type': frag.kind,
                     'length': len(sequence),
                     'sequence_5to3': sequence,
                     'gc_percent': f'{gc_percent(sequence):.1f}',
                     'tm_nn_celsius': f'{melting_temperature(sequence):.1f}',
                     'source_start': frag.source_start,
                     'source_stop': frag.source_stop,
                     'note': frag.note})

    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    with open(fasta_path, 'w') as f:
        for row in rows:
            f.write(f'>{row["name"]} type={row["type"]} length={row["length"]}\n')
            f.write(f'{row["sequence_5to3"]}\n')

    return None


def create_benchling_dna_file(sequence, annotations, file_name, output_name):
    """
    Create a Benchling-readable .dna file.

    :param sequence: str, the full DNA sequence.
    :param annotations: list of tuples, each containing (start, end, label, type) with start and
        end a zero based, half open range of ``sequence``.
    :param file_name: Path, the name of the output .dna file.
    """
    seq_record = SeqRecord(Seq(sequence), id='ExampleID', name=f'{output_name.stem}',
                           description=f'Annotated {output_name}')
    seq_record.annotations['molecule_type'] = 'DNA'
    for start, end, label, feature_type in annotations:
        feature = SeqFeature(FeatureLocation(start, end), type=feature_type,
                             qualifiers={'label': label})
        seq_record.features.append(feature)

    SeqIO.write(seq_record, file_name, 'genbank')


def export_coding_annotation(fragments, args):
    """Writes the coding strand as an annotated Benchling .dna file, 5` -> 3`."""
    coding_fragment = fragments[0]
    sequence = coding_fragment.strand.seq(from_5p=True)
    n_bases = len(sequence)

    # sub_regions are in the strand's 3` -> 5` index space, the .dna file is written 5` -> 3`
    annotations = []
    for start, stop, label in coding_fragment.sub_regions:
        annotations.append((n_bases - stop, n_bases - start, label, 'DNA'))
    annotations.sort()

    dna_path = args.output_name.parent / f'{args.output_name.stem}_coding.dna'
    create_benchling_dna_file(sequence, annotations, dna_path, args.output_name)
    return None


def write_summary(strand, regions, fragments, args):
    """Prints and saves the slicing statistics."""
    staples = [frag for frag in fragments if frag.kind == 'staple']
    staple_lengths = np.array([len(frag.strand) for frag in staples])
    staple_gc = np.array([gc_percent(frag.strand.seq(from_5p=True)) for frag in staples])
    staple_tm = np.array([melting_temperature(frag.strand.seq(from_5p=True)) for frag in staples])

    coding_fragment = fragments[0]
    struct_start, struct_stop = regions['structuring']

    loose = [frag for frag in fragments if frag.kind == 'three_prime']
    loose_note = ''.join(f'\n    Left as its own strand: {frag.name} ({len(frag.strand)} nt)'
                         for frag in loose)

    stat_print = f"""
Structuring strand slicer result info:
    Input structure: {args.topology_file} / {args.strucutre_file}
    Input nucleotides: {len(strand)} in 1 strand
    Output nucleotides: {sum(len(frag.strand) for frag in fragments)} in {len(fragments)} strands

    Region boundaries (oxDNA 3` -> 5` index space):
        3 prime:     {regions['three_prime'][0]} to {regions['three_prime'][1]} ({regions['three_prime'][1] - regions['three_prime'][0]} nt)
        structuring: {struct_start} to {struct_stop} ({struct_stop - struct_start} nt)
        coding:      {regions['coding'][0]} to {regions['coding'][1]} ({regions['coding'][1] - regions['coding'][0]} nt)
        5 prime:     {regions['five_prime'][0]} to {regions['five_prime'][1]} ({regions['five_prime'][1] - regions['five_prime'][0]} nt)

    Coding strand: {len(coding_fragment.strand)} nt{f' ({coding_fragment.note})' if coding_fragment.note else ''}
    Coding/structuring backbone nicked at index {struct_stop - 1} / {struct_stop}{loose_note}

    Slicing method: {args.method}, target staple length {args.staple_length} nt, nicking from the {args.nick_from.replace('_', ' ')}
    Staples: {len(staples)}
    Staple length  min/mean/max: {staple_lengths.min()} / {staple_lengths.mean():.1f} / {staple_lengths.max()} nt
    Staple GC      min/mean/max: {staple_gc.min():.1f} / {staple_gc.mean():.1f} / {staple_gc.max():.1f} %
    Staple Tm (NN) min/mean/max: {np.nanmin(staple_tm):.1f} / {np.nanmean(staple_tm):.1f} / {np.nanmax(staple_tm):.1f} C

    Files written:
""" + '\n'.join(f'        {path}' for path in output_paths(args.output_name)) + '\n'

    print(stat_print)

    with open(args.output_name.with_suffix('.txt'), 'w') as f:
        f.write(stat_print)

    return None


if __name__ == '__main__':
    main()
