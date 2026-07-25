# structuring_strand_slicer

Turns the single-stranded origami written by
[`strand_design/structring_strand_design.py`](../../strand_design/structring_strand_design.py)
into a **coding strand plus a set of DNA origami staples**.

The design script leaves you with one long strand: the 5` half carries the coding sequence (with
optional UTRs), the 3` half is the structuring complement that folds the nanostructure. This
script nicks that strand so the structuring half becomes short, orderable staples.

## What it does

1. **Nicks the coding/structuring backbone bond.** The coding sequence and the 5` UTR come away
   as one strand; the structuring region is left on its own.
2. **Moves the 3` region.** The 3` UTR (stop codon + poly-A) sits at the far end of the
   structuring region in the design output, so nicking would strand it on a staple. It is
   detached and rebuilt onto the coding sequence`s own 3` end — at the coding/structuring nick —
   as a straight helix phased and positioned to continue off that terminus.
3. **Nicks the structuring strand every `N` nucleotides** (default 35), producing the staples.

oxDNA stores a strand 3` -> 5`, so index 0 is the 3` terminus. In that index space the design
output is laid out as:

```
index   0        L3-1  L3                 C0-1  C0           C1-1  C1        N-1
        |  3` UTR  |  |  structuring region  |  |  coding sequence  |  |  5` UTR  |
```

Read 5` -> 3` (descending index) that is: 5` UTR, coding, structuring, 3` UTR.

## Usage

```bash
python3 structuring_strand_slicer.py \
  --strucutre_file  egfp_structring_encoded.dat \
  --topology_file   egfp_structring_encoded.top \
  --coding_sequence egfp_coding_seq.txt \
  --five_prime      fiveprime.txt \
  --three_prime     threeprime.txt \
  --staple_length   35 \
  --output_name     egfp_staples \
  --force_overwrite
```

[`run_structuring_strand_slicer.sh`](run_structuring_strand_slicer.sh) runs exactly this against
the [example/](../../example) design.

The three sequence files must be the ones used for the design run — they are what locates the
region boundaries, and the script verifies the structure`s bases against them before slicing. If
you no longer have them, pass `--coding_length` / `--five_prime_length` / `--three_prime_length`
instead.

### Arguments

| Flag | Required | Description |
|------|:---:|-------------|
| `-s`, `--strucutre_file` | yes | oxDNA `.dat` written by the design script |
| `-t`, `--topology_file` | yes | oxDNA `.top` written by the design script |
| `-c`, `--coding_sequence` | yes* | `.txt`/`.fasta` coding sequence used for the design |
| `--five_prime`, `--three_prime` | no | `.txt`/`.fasta` UTRs used for the design |
| `--coding_length`, `--five_prime_length`, `--three_prime_length` | no | region lengths, instead of the files above |
| `-n`, `--staple_length` | no | nucleotides between nicks (default 35) |
| `--min_staple_length` | no | shortest staple allowed; a shorter leftover is merged into its neighbour (default 15) |
| `-m`, `--method` | no | `uniform` (cut exactly every `N`) or `even` (near-equal staples) |
| `--nick_from` | no | which end the fixed-length cuts start from: `coding_junction` (default) or `three_prime_junction` |
| `--three_prime_placement` | no | `coding_3p` (default, move it) or `keep` (leave it as its own strand) |
| `--coding_as_rna` | no | report the coding strand as RNA (T -> U) in the sequence outputs |
| `--no_verify_regions` | no | skip the base-level check against the sequence files |
| `-o`, `--output_name` | no | output basename (default `structuring_staples`) |
| `--force_overwrite` | no | overwrite existing output files |

\* or `--coding_length`.

### Outputs

For an `--output_name` of `NAME`:

| File | Contents |
|------|----------|
| `NAME.top`, `NAME.dat` | Sliced structure in oxDNA format |
| `NAME.oxview` | oxView file, coloured per staple and per region, with a selections legend |
| `NAME_strands.csv` | Primary sequences 5` -> 3` with length, GC%, nearest-neighbour Tm, and the source index range |
| `NAME_strands.fasta` | The same sequences as FASTA |
| `NAME_coding.dna` | Benchling-readable (GenBank) coding strand annotated with 5` UTR / coding / 3` UTR |
| `NAME.txt` | Summary statistics |

## Adding another slicing method

`chunk_lengths()` turns the structuring region length into a list of staple lengths, and
`structuring_bounds()` turns those into index ranges. A new strategy — nicking at duplex domain
boundaries, for instance, using the `hb_list_traj.txt` the design script writes — only needs a
new branch in `chunk_lengths()` (or a replacement for `structuring_bounds()`) plus an entry in
`SLICING_METHODS`.

## Dependencies

`ipy_oxdna`, `numpy`, `matplotlib` and `biopython`. Unlike the design script this one does not
need `oxpy`, `nupack` or the `oat` CLI.
