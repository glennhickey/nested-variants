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
# 4. pantree-extract.py on yeast test VCF
############################################################################
PANTREE_VCF="$SCRIPT_DIR/../yeast-test/pantree-chrI.vcf"
if command -v bcftools >/dev/null 2>&1 && command -v python3 >/dev/null 2>&1 && [ -f "$PANTREE_VCF" ]; then
    echo ""
    echo "=== pantree-extract.py test ==="
    PTDIR=$(mktemp -d)
    if python3 "$SCRIPT_DIR/pantree-extract.py" --vcf "$PANTREE_VCF" --output "$PTDIR/pantree.records.tsv" 2>&1; then
        # Verify output has chrI chromosomes (not chr1)
        if head -2 "$PTDIR/pantree.records.tsv" | grep -q "^chrI"; then
            echo "  PASS  pantree-extract.py (chrI records)"
        else
            echo "  FAIL  pantree-extract.py (expected chrI in output)"
            FAIL=1
        fi
        # Verify no chr1 records leaked through
        if grep -q "^chr1[^I]" "$PTDIR/pantree.records.tsv" 2>/dev/null; then
            echo "  FAIL  pantree-extract.py (unexpected chr1 records)"
            FAIL=1
        fi
    else
        echo "  FAIL  pantree-extract.py (non-zero exit)"
        FAIL=1
    fi
    rm -rf "$PTDIR"
fi

############################################################################
echo ""
if [ "$FAIL" -eq 0 ]; then
    echo "All tests passed."
else
    echo "Some tests FAILED."
    exit 1
fi
