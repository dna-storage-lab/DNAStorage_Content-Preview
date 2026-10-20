#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
PROJECT_DIR=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
BUILD_DIR="$SCRIPT_DIR/build"

ALIGNED_ROOT=${1:-"$PROJECT_DIR/result/aligned_reads"}
OUTPUT_ROOT=${2:-"$PROJECT_DIR/result/fb_results"}
THREAD_COUNT=${3:-4}

mkdir -p "$BUILD_DIR" "$OUTPUT_ROOT"

if [[ ! -f "$SCRIPT_DIR/HalfMaker.h" ]]; then
    echo "Missing $SCRIPT_DIR/HalfMaker.h; it is required to compile Step 3." >&2
    exit 1
fi

echo "Compiling Step 3 (forward-backward inference)..."
g++ -std=c++17 -O2 -Wall -Wextra -pedantic -pthread \
    "$SCRIPT_DIR/fb_probability_extraction.cpp" \
    "$SCRIPT_DIR/fb_indel_correction.cpp" \
    -o "$BUILD_DIR/fb_indel_correction"

while IFS= read -r -d '' match_file; do
    relative_path=${match_file#"$ALIGNED_ROOT/"}
    relative_dir=$(dirname -- "$relative_path")
    file_name=$(basename -- "$relative_path")
    sample_name=${file_name%_match_reads_lower.txt}
    output_dir="$OUTPUT_ROOT/$relative_dir"
    output_prefix="$output_dir/$sample_name"

    mkdir -p "$output_dir"
    "$BUILD_DIR/fb_indel_correction" \
        "$match_file" \
        "$output_prefix" \
        "$THREAD_COUNT"
done < <(find "$ALIGNED_ROOT" -type f -name '*_match_reads_lower.txt' -print0 | sort -z)

echo "Step 3 complete: $OUTPUT_ROOT"
