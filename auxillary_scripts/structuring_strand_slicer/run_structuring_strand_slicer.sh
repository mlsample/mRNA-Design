#!/usr/bin/env bash
# Slices the example design (../../example) into DNA origami staples.
# Run ../../example/run_strand_design.sh first, it writes the .top/.dat this reads.
set -euo pipefail

example_dir="$(dirname "$0")/../../example"

python3 "$(dirname "$0")/structuring_strand_slicer.py" \
 --strucutre_file "${example_dir}/egfp_structring_encoded.dat" \
 --topology_file  "${example_dir}/egfp_structring_encoded.top" \
 --coding_sequence "${example_dir}/egfp_coding_seq.txt" \
 --five_prime  "${example_dir}/fiveprime.txt" \
 --three_prime "${example_dir}/threeprime.txt" \
 --staple_length 35 \
 --output_name egfp_staples \
 --force_overwrite
