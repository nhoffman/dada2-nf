#!/usr/bin/env bash
# Run nextflow pipeline for each test params.json file.
# Usage: ./test/run_tests.sh [test_name ...]
# If no args given, runs all tests found in test/**/params.json.
# Exit status is the number of failed tests.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

PASS=0
FAIL=0
FAILURES=()

run_test() {
    local params_file="$1"
    local test_name
    test_name="$(realpath --relative-to="$SCRIPT_DIR" "$params_file" | sed 's|/params.json||')"

    echo "========================================"
    echo "TEST: $test_name"
    echo "  params: $params_file"
    echo "========================================"

    if nextflow run "$PROJECT_DIR/main.nf" \
            -params-file "$params_file" \
            -ansi-log false \
            2>&1; then
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

if [[ $# -gt 0 ]]; then
    # Run only the specified tests
    for name in "$@"; do
        params="$SCRIPT_DIR/$name/params.json"
        if [[ ! -f "$params" ]]; then
            echo "ERROR: params.json not found for test '$name' at $params" >&2
            ((FAIL++)) || true
            FAILURES+=("$name")
        else
            run_test "$params"
        fi
    done
else
    # Run all tests
    while IFS= read -r params_file; do
        run_test "$params_file"
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
