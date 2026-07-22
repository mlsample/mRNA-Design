# nanostructured_expressing_nucleic_acids

Design the sequence of a **single-stranded origami (ssOrigami)** so that it both folds into a
target nanostructure *and* carries a translatable protein-coding sequence.

A single-stranded RNA (or DNA) origami folds into a compact nanostructure through internal,
intramolecular base pairing. This package takes an already-folded ssOrigami — defined by an
oxDNA structure and its measured base-pairing pattern — and rewrites its sequence so that a
chosen coding sequence (e.g. eGFP) is embedded in one half of the strand while the other half
is re-mutated to remain the Watson–Crick complement of the coding region. The result is a
strand whose fold is preserved but whose sequence now encodes a protein.

The core algorithm is the **structuring-strand design** in
[strand_design/structring_strand_design.py](strand_design/structring_strand_design.py).

## How it works

Given a folded ssOrigami, the tool:

1. **Reads the fold.** Runs `oat output_bonds` (from `oxDNA_analysis_tools`) over the structure
   (or a trajectory) to recover which nucleotides are hydrogen-bonded to which, building a
   `pair_map` of the actual base pairs that hold the nanostructure together.
2. **Places the coding sequence.** Splits the strand in half and lays the coding sequence into
   the second half (in reverse index order, so the structure's 3'→5' direction reads 5'→3' as
   coding), then mutates those nucleotides to the requested sequence.
3. **Repairs the complement.** For every structuring nucleotide paired to a coding nucleotide,
   it mutates the partner to the Watson–Crick complement so the duplex — and therefore the fold —
   is retained. It handles the awkward cases explicitly: coding that pairs with coding, coding
   with no pairing partner, and unpaired structuring regions.
4. **Validates.** Asserts that the embedded coding sequence matches the request and that every
   preserved pair is still complementary.
5. **Scores destabilization.** Groups the pairs into duplex domains and compares each domain's
   NUPACK free energy before and after mutation, reporting per-domain percent energy change and
   flagging strongly destabilized domains.
6. **Adds UTRs (optional).** Prepends a 3' sequence and appends a 5' sequence (e.g. Kozak/start
   and stop/poly-A) to the strand.
7. **Exports** the redesigned structure in several formats (see below).

### Outputs

For an `--output_name` of `NAME`, the tool writes:

| File | Contents |
|------|----------|
| `NAME.dat`, `NAME.top` | Redesigned structure in oxDNA format |
| `NAME.oxview` | oxView file, colored by region with a selections legend (coding, structuring, destabilized domains, UTRs, …) |
| `NAME.dna` | Benchling-readable (GenBank) file with region annotations |
| `NAME.png` | Histogram of per-domain percent free-energy change |
| `NAME.txt` | Summary statistics (pairing %, domain count, mean energy change, …) |

## Installation

```bash
pip install -e .
```

This installs the `strand_design` package (`numpy`, `pandas`, `matplotlib`, `seaborn`,
`biopython`). It additionally depends on tools that are **not** installed automatically and must
be available in your environment:

- [`oxDNA_analysis_tools`](https://github.com/lorenzo-rovigatti/oxDNA) (the `oat` CLI) and `oxpy`
- `ipy_oxdna` (structure editing / oxView export)
- [`nupack`](https://nupack.org) (domain free-energy scoring, used by
  [strand_design/RNA_NN_better.py](strand_design/RNA_NN_better.py))

## Usage

```bash
python3 strand_design/structring_strand_design.py \
  --strucutre_file  ssOrigami_last_conf.dat \
  --topology_file   ssOrigami.top \
  --input_md_file   inputMD \
  --traj_file       trajectory.dat \
  --coding_sequence coding_seq.txt \
  --five_prime      fiveprime.txt \
  --three_prime     threeprime.txt \
  --output_name     my_design \
  --force_overwrite
```

### Arguments

| Flag | Required | Description |
|------|:---:|-------------|
| `-s`, `--strucutre_file` | yes | oxDNA `.dat` configuration of the folded ssOrigami |
| `-t`, `--topology_file` | yes | oxDNA `.top` topology file |
| `-i`, `--input_md_file` | yes | oxDNA input file, used by `output_bonds` to compute base pairs |
| `-c`, `--coding_sequence` | yes | `.txt` (one line) or `.fasta` file of the coding sequence, 5'→3' (A/T/C/G/U) |
| `--traj_file` | no | Trajectory `.dat`; when given, base pairs are averaged over the trajectory |
| `--five_prime` | no | `.txt`/`.fasta` of bases (5'→3') to append to the 5' end (e.g. Kozak + start codon) |
| `--three_prime` | no | `.txt`/`.fasta` of bases (5'→3') to prepend to the 3' end (e.g. stop codon + poly-A) |
| `-o`, `--output_name` | no | Output basename (default: `coding_sequence_embbeded_strucutre`) |
| `--force_overwrite` | no | Overwrite existing output files |

The coding sequence must be no longer than half of the structure.

## Example

The [example/](example/) directory contains a complete, runnable case: an eGFP coding sequence
embedded into a `science`-style ssRNA origami, with 5'/3' UTRs. Run it with:

```bash
cd example
bash run_strand_design.sh
```

See [example/run_strand_design.sh](example/run_strand_design.sh) for the exact invocation.

## Repository layout

| Path | Description |
|------|-------------|
| [strand_design/structring_strand_design.py](strand_design/structring_strand_design.py) | Main CLI: embed a coding sequence into an ssOrigami |
| [strand_design/folder_base.py](strand_design/folder_base.py) | `ssOrigamiParse` class; base-pair stretch decoding and per-domain energy analysis |
| [strand_design/RNA_NN_better.py](strand_design/RNA_NN_better.py) | NUPACK free-energy helper for duplex domains |
| [example/](example/) | Runnable eGFP example with input files and a driver script |
| [auxillary_scripts/](auxillary_scripts/) | JavaScript ports and notebooks for mutating/exporting structures |
| [tests/](tests/) | `pytest` suite for the design pipeline |

## Tests

```bash
pytest
```

## License

MIT — see [LICENSE](LICENSE).
