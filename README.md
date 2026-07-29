# mRNA-Design

Design the sequence of a **single-stranded origami (ssOrigami)** so that it both folds into a
target nanostructure *and* carries a translatable protein-coding sequence.

## Background

A single-stranded RNA (or DNA) origami folds into a compact nanostructure through internal,
intramolecular base pairing.

mRNA-Design starts from an already-folded ssOrigami — defined by an oxDNA structure and its
measured base-pairing pattern — and rewrites its sequence so that:

- a chosen coding sequence (e.g. eGFP) is embedded in one half of the strand, and
- the other half is re-mutated to remain the Watson–Crick complement of the coding region.

The result is a strand whose fold is preserved but whose sequence now encodes a protein. The
core algorithm is the **structuring-strand design** in
[strand_design/structring_strand_design.py](strand_design/structring_strand_design.py).

## Contents

- [How it works](#how-it-works)
- [Installation](#installation)
- [Usage](#usage)
- [Example](#example)
- [Repository layout](#repository-layout)
- [Tests](#tests)

## How it works

Given a folded ssOrigami, the tool:

1. **Reads the fold.** Runs `oat output_bonds` (from `oxDNA_analysis_tools`) over the structure —
   or over a trajectory, in which case hydrogen-bond energies are averaged across frames. Pairs
   with an HB energy below `-0.1` are kept, and if a nucleotide ends up with more than one partner
   only the strongest pair survives. The result is a `pair_map` of the base pairs that actually
   hold the nanostructure together.
2. **Places the coding sequence.** Splits the strand in half and lays the coding sequence into the
   upper-index half, descending from the last index. oxDNA stores the strand 3'→5', so this puts
   the coding sequence at the 5' end reading 5'→3'. Those nucleotides are then mutated to the
   requested sequence.
3. **Repairs the complement.** For every structuring nucleotide paired to a coding nucleotide,
   it mutates the partner to the Watson–Crick complement so the duplex — and therefore the fold —
   is retained. It handles the awkward cases explicitly: coding that pairs with coding, coding
   with no pairing partner, and unpaired structuring regions.
4. **Validates.** Asserts that the embedded coding sequence matches the request and that every
   preserved pair is still complementary.
5. **Scores destabilization.** Groups the pairs into duplex domains — runs of at least 4
   consecutive base pairs — and compares each domain's NUPACK free energy (RNA model, 37 °C)
   before and after mutation. Domains that lose more than 40% of their free energy are flagged as
   strongly destabilized.
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

It also leaves the extracted base-pair list in the working directory as `hblist.txt` (raw
`output_bonds` output) and `hb_list_traj.csv` / `hb_list_traj.txt` (the filtered, trajectory-averaged
pairs). These are inputs to the energy analysis, and are useful for inspecting the fold directly.

## Installation

```bash
pip install -e .
```

This pulls in the PyPI dependencies: `numpy`, `pandas`, `matplotlib`, `seaborn`, and `biopython`.

### Installed separately

The tool also needs the following, which are not on PyPI and are **not** installed by the command
above:

| Dependency | Used for |
|------------|----------|
| [`oxDNA_analysis_tools`](https://github.com/lorenzo-rovigatti/oxDNA) (the `oat` CLI) and `oxpy` | Extracting base pairs from the structure |
| `ipy_oxdna` | Structure editing / oxView export |
| [`nupack`](https://nupack.org) | Domain free-energy scoring, via [strand_design/RNA_NN_better.py](strand_design/RNA_NN_better.py) |
| `gawk` | Filtering the `output_bonds` output; called as a shell command |

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

Three things to watch for:

- **Run it from the directory holding your input files.** `output_bonds` writes its `hblist.txt`
  next to the oxDNA input file, but the script reads it back from the current working directory, so
  the two must be the same place. This is why the example below does `cd example` first.
- **The coding sequence must be no longer than half of the structure.**
- **Sequence files must be a single line with no trailing newline.** A stray newline is read as
  part of the sequence and trips the A/T/C/G/U validation. This applies to `--coding_sequence`,
  `--five_prime`, and `--three_prime`; for `.fasta` input, exactly two lines are expected.

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

## Example

The [example/](example/) directory contains a complete, runnable case: an eGFP coding sequence
embedded into a `science`-style ssRNA origami, with 5'/3' UTRs. Run it with:

```bash
cd example
bash run_strand_design.sh
```

See [example/run_strand_design.sh](example/run_strand_design.sh) for the exact invocation. The base
pairs come from `trajectory.dat` rather than a single configuration, and `inputMD` is an `RNA2`
sequence-dependent input, so `rna_sequence_dependent_parameters.txt` must sit alongside it.

## Repository layout

| Path | Description |
|------|-------------|
| [strand_design/structring_strand_design.py](strand_design/structring_strand_design.py) | Main CLI: embed a coding sequence into an ssOrigami |
| [strand_design/folder_base.py](strand_design/folder_base.py) | `ssOrigamiParse` class; base-pair stretch decoding and per-domain energy analysis |
| [strand_design/RNA_NN_better.py](strand_design/RNA_NN_better.py) | NUPACK free-energy helper for duplex domains |
| [example/](example/) | Runnable eGFP example with input files and a driver script |
| [Example_ssRNA_Origamis/](Example_ssRNA_Origamis/) | Additional folded ssRNA origami structures (ACS Nano square, Science rectangles) to design against |
| [tests/](tests/) | `pytest` suite for the design pipeline |

## Tests

```bash
pytest
```

Note that the test fixtures read structures from a `strucutres/` directory that is not part of the
repository, so the suite does not run from a fresh clone as-is.

## License

GPL-3.0 — see [LICENSE](LICENSE).
