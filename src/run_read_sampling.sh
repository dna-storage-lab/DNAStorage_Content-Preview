#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
PROJECT_DIR=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
BUILD_DIR="$SCRIPT_DIR/build"

FASTQ_DIR=${1:-"$PROJECT_DIR/sequencing_data"}
OUTPUT_ROOT=${2:-"$PROJECT_DIR/result/sampled_fastq"}
REPLICATE_COUNT=${3:-100}
START_READ=${4:-1}

mkdir -p "$BUILD_DIR" "$OUTPUT_ROOT"

echo "Compiling Step 1 (FASTQ read sampling)..."
g++ -std=c++11 -O2 -Wall -Wextra -pedantic \
    "$SCRIPT_DIR/read_sampling.cpp" \
    -o "$BUILD_DIR/read_sampling"

sample_dataset() {
    local dataset=$1
    local fastq_name=$2
    local max_reads=$3
    local input_fastq="$FASTQ_DIR/$fastq_name"

    if [[ ! -f "$input_fastq" ]]; then
        echo "Missing FASTQ input: $input_fastq" >&2
        return 1
    fi

    echo "Sampling $dataset from $input_fastq (1..$max_reads reads per replicate)"
    local reads_per_replicate
    for ((reads_per_replicate = 1; reads_per_replicate <= max_reads; ++reads_per_replicate)); do
        local output_dir="$OUTPUT_ROOT/$dataset/coverage_$reads_per_replicate"
        mkdir -p "$output_dir"
        "$BUILD_DIR/read_sampling" \
            "$input_fastq" \
            "$output_dir" \
            "$reads_per_replicate" \
            "$REPLICATE_COUNT" \
            "$START_READ"
    done
}

sample_dataset poem1   poem1.fastq   25
# sample_dataset poem2   poem2.fastq   25
# sample_dataset poem3   poem3.fastq   25
# sample_dataset poem4_5 'poem4&5.fastq' 25
# sample_dataset poem6   poem6.fastq   50
# sample_dataset poem7   poem7.fastq   25
# sample_dataset poem8   poem8.fastq   25

echo "Step 1 complete: $OUTPUT_ROOT"
