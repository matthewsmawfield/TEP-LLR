#!/bin/bash
# CI Check: Verify manuscript values match pipeline outputs before PDF generation
# Usage: ./scripts/utils/ci_check.sh

set -e

echo "=========================================="
echo "TEP-LLR Manuscript Integrity CI Check"
echo "=========================================="
echo ""

# Run verification script
python scripts/utils/verify_value_consistency.py
EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ CI Check PASSED - Manuscript values match pipeline outputs"
    echo "  PDF generation can proceed safely."
    exit 0
else
    echo "❌ CI Check FAILED - Discrepancies detected between manuscript and pipeline"
    echo "  Fix the issues above before generating PDF."
    exit 1
fi
