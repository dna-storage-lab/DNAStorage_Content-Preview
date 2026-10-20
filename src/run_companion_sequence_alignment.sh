#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
PROJECT_DIR=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
BUILD_DIR="$SCRIPT_DIR/build"

SAMPLED_ROOT=${1:-"$PROJECT_DIR/result/sampled_fastq"}
COMPANION_SEQUENCE_DIR=${2:-"$PROJECT_DIR/companion_sequence"}
OUTPUT_ROOT=${3:-"$PROJECT_DIR/result/aligned_reads"}

mkdir -p "$BUILD_DIR" "$OUTPUT_ROOT"

if [[ ! -f "$SCRIPT_DIR/edlib.h" ]]; then
    echo "Missing $SCRIPT_DIR/edlib.h. Install/copy the Edlib header before compiling." >&2
    exit 1
fi

echo "Compiling Step 2 (demapping and companion alignment)..."
g++ -std=c++11 -O2 -Wall -Wextra -pedantic \
    "$SCRIPT_DIR/edlib.cpp" \
    "$SCRIPT_DIR/companion_sequence_alignment.cpp" \
    -o "$BUILD_DIR/companion_sequence_alignment"

align_target() {
    local target=$1
    local sampled_dataset=$2
    local companion_sequence_name=$3
    local sampled_dir="$SAMPLED_ROOT/$sampled_dataset"
    local companion_sequence_file="$COMPANION_SEQUENCE_DIR/$companion_sequence_name"

    if [[ ! -d "$sampled_dir" ]]; then
        echo "Missing sampled FASTQ directory: $sampled_dir" >&2
        return 1
    fi
    if [[ ! -f "$companion_sequence_file" ]]; then
        echo "Missing companion sequence: $companion_sequence_file" >&2
        return 1
    fi

    echo "Aligning $target using $companion_sequence_name"
    while IFS= read -r -d '' sample_file; do
        local relative_path=${sample_file#"$sampled_dir/"}
        local relative_dir
        relative_dir=$(dirname -- "$relative_path")
        local sample_name
        sample_name=$(basename -- "$relative_path" .fastq)
        local target_dir="$OUTPUT_ROOT/$target/$relative_dir"
        local match_file="$target_dir/${sample_name}_match_reads_lower.txt"

        mkdir -p "$target_dir"
        "$BUILD_DIR/companion_sequence_alignment" "$companion_sequence_file" "$sample_file" "$match_file"
    done < <(find "$sampled_dir" -type f -name '*.fastq' -print0 | sort -z)
}

align_target poem1 poem1   companion_sequence_poem1.txt
# align_target poem2 poem2   companion_sequence_poem2.txt
# align_target poem3 poem3   companion_sequence_poem3.txt
# align_target poem4 poem4_5 companion_sequence_poem4.txt
# align_target poem5 poem4_5 companion_sequence_poem5.txt
# align_target poem6 poem6   companion_sequence_poem6.txt
# align_target poem7 poem7   companion_sequence_poem7.txt
# align_target poem8 poem8   companion_sequence_poem8.txt

echo "Step 2 complete: $OUTPUT_ROOT"
