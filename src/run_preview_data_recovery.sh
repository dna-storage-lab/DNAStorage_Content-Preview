#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
PROJECT_DIR=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
BUILD_DIR="$SCRIPT_DIR/build"

FB_ROOT=${1:-"$PROJECT_DIR/result/fb_results"}
OUTPUT_ROOT=${2:-"$PROJECT_DIR/result/recovered"}

mkdir -p "$BUILD_DIR" "$OUTPUT_ROOT"

echo "Compiling Step 4 (probability merging and data recovery)..."
g++ -std=c++11 -O2 -Wall -Wextra -pedantic \
    "$SCRIPT_DIR/probability_merging.cpp" \
    -o "$BUILD_DIR/probability_merging"
g++ -std=c++11 -O2 -Wall -Wextra -pedantic \
    "$SCRIPT_DIR/preview_data_recovery.cpp" \
    -o "$BUILD_DIR/preview_data_recovery"

recover_target() {
    local target=$1
    local bit_length=$2
    local target_root="$FB_ROOT/$target"

    if [[ ! -d "$target_root" ]]; then
        echo "Missing Step 3 directory: $target_root" >&2
        return 1
    fi

    echo "Recovering $target ($bit_length information bits)"
    while IFS= read -r -d '' prob0_file; do
        local relative_path=${prob0_file#"$target_root/"}
        local relative_dir
        relative_dir=$(dirname -- "$relative_path")
        local file_name
        file_name=$(basename -- "$relative_path")
        local sample_name=${file_name%_prob0.txt}
        local prob1_file="${prob0_file%_prob0.txt}_prob1.txt"
        local output_dir="$OUTPUT_ROOT/$target/$relative_dir"
        local voting_result="$output_dir/${sample_name}_voting.txt"
        local recovered_file="$output_dir/${sample_name}_preview.txt"

        if [[ ! -s "$prob0_file" || ! -s "$prob1_file" ]]; then
            echo "Skipping $target/$relative_path: no aligned-read probabilities." >&2
            continue
        fi

        mkdir -p "$output_dir"
        "$BUILD_DIR/probability_merging" \
            "$prob0_file" \
            "$prob1_file" \
            "$voting_result" \
            "$bit_length"
        "$BUILD_DIR/preview_data_recovery" \
            "$voting_result" \
            "$recovered_file" \
            "$bit_length"
    done < <(find "$target_root" -type f -name '*_prob0.txt' -print0 | sort -z)
}

recover_target poem1 1000
# recover_target poem2 2016
# recover_target poem3 360
# recover_target poem4 992
# recover_target poem5 264
# recover_target poem6 1000
# recover_target poem7 816
# recover_target poem8 1224

echo "Step 4 complete: $OUTPUT_ROOT"
