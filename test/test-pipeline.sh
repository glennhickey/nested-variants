#!/bin/bash
# test-pipeline.sh — CI test script for nested-variants
# Runs shellcheck on all shell scripts and verifies --help exits cleanly.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")/.." && pwd)/scripts"
FAIL=0

############################################################################
# 1. shellcheck all .sh scripts
############################################################################
echo "=== shellcheck ==="
if command -v shellcheck >/dev/null 2>&1; then
    for script in "$SCRIPT_DIR"/*.sh; do
        name=$(basename "$script")
        if shellcheck -x "$script"; then
            echo "  PASS  $name"
        else
            echo "  FAIL  $name"
            FAIL=1
        fi
    done
else
    echo "  SKIP  shellcheck not installed"
fi

############################################################################
# 2. --help exits 0 for every shell script
############################################################################
echo ""
echo "=== --help flag ==="
for script in "$SCRIPT_DIR"/*.sh; do
    name=$(basename "$script")
    if bash "$script" --help >/dev/null 2>&1; then
        echo "  PASS  $name --help"
    else
        echo "  FAIL  $name --help (exit $?)"
        FAIL=1
    fi
done

############################################################################
# 3. Optional: split-ref.sh with a tiny VCF fixture
############################################################################
if command -v bcftools >/dev/null 2>&1 && [ -f "$SCRIPT_DIR/../data/test.vcf.gz" ]; then
    echo ""
    echo "=== split-ref.sh integration test ==="
    TMPDIR=$(mktemp -d)
    trap 'rm -rf "$TMPDIR"' EXIT
    cd "$TMPDIR"
    if bash "$SCRIPT_DIR/split-ref.sh" -v "$SCRIPT_DIR/../data/test.vcf.gz" -p GRCh38; then
        echo "  PASS  split-ref.sh integration"
    else
        echo "  FAIL  split-ref.sh integration"
        FAIL=1
    fi
fi

############################################################################
echo ""
if [ "$FAIL" -eq 0 ]; then
    echo "All tests passed."
else
    echo "Some tests FAILED."
    exit 1
fi
