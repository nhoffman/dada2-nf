#!/usr/bin/env bash
# Run nextflow pipeline for each test params.json file.
# Usage: ./test/run_tests.sh [test_name ...] [-- nextflow_arg ...]
# If no args given, runs all tests found in test/**/params.json.
# Exit status is the number of failed tests.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

PASS=0
FAIL=0
FAILURES=()

params_output_dir() {
    python3 - "$1" <<'PY'
import json
import sys

with open(sys.argv[1]) as handle:
    params = json.load(handle)
print(params.get('output', 'output'))
PY
}

check_base_files() {
    local output_dir="$1"
    local test_dir="$2"
    local base_file="$test_dir/base-files.sha256"

    if [[ ! -f "$base_file" ]]; then
        echo "ERROR: baseline file not found: $base_file" >&2
        return 1
    fi

    if [[ ! -s "$base_file" ]]; then
        echo "No baseline checks configured: $base_file"
        return 0
    fi

    (cd "$output_dir" && sha256sum -c "$base_file")
}

run_test() {
    local params_file="$1"
    shift
    local test_name
    local test_dir
    local output_dir
    test_name="$(realpath --relative-to="$SCRIPT_DIR" "$params_file" | sed 's|/params.json||')"
    test_dir="$(dirname "$params_file")"
    output_dir="$(params_output_dir "$params_file")"

    echo "========================================"
    echo "TEST: $test_name"
    echo "  params: $params_file"
    echo "  output: $output_dir"
    echo "========================================"

    if nextflow run "$PROJECT_DIR/main.nf" \
            -params-file "$params_file" \
            -ansi-log false \
            "$@" \
            2>&1 && check_base_files "$output_dir" "$test_dir"; then
        echo "PASSED: $test_name"
        ((PASS++)) || true
    else
        echo "FAILED: $test_name"
        ((FAIL++)) || true
        FAILURES+=("$test_name")
    fi
    echo
}

cd "$PROJECT_DIR"

TEST_ARGS=()
NEXTFLOW_ARGS=()
while [[ $# -gt 0 ]]; do
    if [[ "$1" == "--" ]]; then
        shift
        NEXTFLOW_ARGS=("$@")
        break
    fi
    TEST_ARGS+=("$1")
    shift
done

if [[ ${#TEST_ARGS[@]} -gt 0 ]]; then
    # Run only the specified tests
    for name in "${TEST_ARGS[@]}"; do
        params="$SCRIPT_DIR/$name/params.json"
        if [[ ! -f "$params" ]]; then
            echo "ERROR: params.json not found for test '$name' at $params" >&2
            ((FAIL++)) || true
            FAILURES+=("$name")
        else
            run_test "$params" "${NEXTFLOW_ARGS[@]}"
        fi
    done
else
    # Run all tests
    while IFS= read -r params_file; do
        run_test "$params_file" "${NEXTFLOW_ARGS[@]}"
    done < <(find "$SCRIPT_DIR" -name "params.json" | sort)
fi

echo "========================================"
echo "Results: $PASS passed, $FAIL failed"
if [[ ${#FAILURES[@]} -gt 0 ]]; then
    echo "Failed tests:"
    for f in "${FAILURES[@]}"; do
        echo "  - $f"
    done
fi
echo "========================================"

exit "$FAIL"
