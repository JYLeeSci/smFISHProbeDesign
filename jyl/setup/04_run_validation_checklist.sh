#!/bin/bash
# Validation Checklist Runner
# Purpose: Run validation checklist items from jyl_implementation_plan.md Section 6.5
# Output: Logs results to jyl/setup/logs/validation_checklist.log

set -euo pipefail

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Paths
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LOG_FILE="$REPO_ROOT/jyl/setup/logs/validation_checklist.log"

# Initialize log
echo "=== Validation Checklist ===" > "$LOG_FILE"
echo "Date: $(date)" >> "$LOG_FILE"
echo "Location: $REPO_ROOT" >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

# Results tracking
PASS_COUNT=0
FAIL_COUNT=0
TOTAL_COUNT=5

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Validation Checklist (Section 6.5)${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Check 1: probedesign --version
echo -e "${BLUE}[1/5] Checking probedesign version...${NC}"
echo "[1/5] probedesign --version" >> "$LOG_FILE"
VERSION_OUTPUT=$(probedesign --version 2>&1 || echo "ERROR")
if echo "$VERSION_OUTPUT" | grep -q "0.1.0"; then
    echo -e "${GREEN}✓ PASS: probedesign version 0.1.0${NC}"
    echo "Result: PASS - $VERSION_OUTPUT" >> "$LOG_FILE"
    ((PASS_COUNT++))
else
    echo -e "${RED}✗ FAIL: Expected 0.1.0, got: $VERSION_OUTPUT${NC}"
    echo "Result: FAIL - $VERSION_OUTPUT" >> "$LOG_FILE"
    ((FAIL_COUNT++))
fi
echo "" >> "$LOG_FILE"
echo ""

# Check 2: bowtie --version
echo -e "${BLUE}[2/5] Checking bowtie version...${NC}"
echo "[2/5] bowtie --version" >> "$LOG_FILE"
BOWTIE_OUTPUT=$(bowtie --version 2>&1 | head -1 || echo "ERROR")
if echo "$BOWTIE_OUTPUT" | grep -q "bowtie-align-s version 1.3.1"; then
    echo -e "${GREEN}✓ PASS: bowtie-align-s version 1.3.1${NC}"
    echo "Result: PASS - $BOWTIE_OUTPUT" >> "$LOG_FILE"
    ((PASS_COUNT++))
else
    echo -e "${RED}✗ FAIL: Expected bowtie-align-s 1.3.1, got: $BOWTIE_OUTPUT${NC}"
    echo "Result: FAIL - $BOWTIE_OUTPUT" >> "$LOG_FILE"
    ((FAIL_COUNT++))
fi
echo "" >> "$LOG_FILE"
echo ""

# Check 3: CDKN1A_32 test
echo -e "${BLUE}[3/5] Running CDKN1A_32 test (32 probes expected)...${NC}"
echo "[3/5] CDKN1A_32 test" >> "$LOG_FILE"
cd "$REPO_ROOT"

OUTPUT_DIR="/tmp/validation_cdkn1a"
rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

probedesign design test_cases/CDKN1A_32/CDKN1A.fa \
  --probes 32 \
  --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa \
  -o "$OUTPUT_DIR/test" > /dev/null 2>> "$LOG_FILE"

ACTUAL_COUNT=$(wc -l < "$OUTPUT_DIR/test_oligos.txt" | tr -d ' ')
if [ "$ACTUAL_COUNT" -eq 32 ]; then
    # Compare sequences
    awk -F'\t' '{print $(NF-1)}' "$OUTPUT_DIR/test_oligos.txt" | sort > /tmp/val_actual.txt
    awk -F'\t' '{print $(NF-1)}' test_cases/CDKN1A_32/CDKN1A_32_genomemaskoff_oligos.txt | sort > /tmp/val_expected.txt

    if diff -q /tmp/val_actual.txt /tmp/val_expected.txt > /dev/null 2>&1; then
        echo -e "${GREEN}✓ PASS: 32/32 probes match (100%)${NC}"
        echo "Result: PASS - 32/32 probes match" >> "$LOG_FILE"
        ((PASS_COUNT++))
    else
        MATCHING=$(comm -12 /tmp/val_actual.txt /tmp/val_expected.txt | wc -l | tr -d ' ')
        echo -e "${RED}✗ FAIL: Only $MATCHING/32 probes match${NC}"
        echo "Result: FAIL - $MATCHING/32 probes match" >> "$LOG_FILE"
        ((FAIL_COUNT++))
    fi
else
    echo -e "${RED}✗ FAIL: Expected 32 probes, got $ACTUAL_COUNT${NC}"
    echo "Result: FAIL - Expected 32 probes, got $ACTUAL_COUNT" >> "$LOG_FILE"
    ((FAIL_COUNT++))
fi
echo "" >> "$LOG_FILE"
echo ""

# Check 4: KRT19 test (with genome indices)
echo -e "${BLUE}[4/5] Running KRT19 test (6 probes expected, requires indices)...${NC}"
echo "[4/5] KRT19 test" >> "$LOG_FILE"

OUTPUT_DIR="/tmp/validation_krt19"
rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

if probedesign design test_cases/KRT19_withUTRs/KRT19_withUTRs.fa \
  --probes 32 \
  --pseudogene-mask \
  --genome-mask \
  --index-dir bowtie_indexes \
  -o "$OUTPUT_DIR/test" > /dev/null 2>> "$LOG_FILE"; then

    ACTUAL_COUNT=$(wc -l < "$OUTPUT_DIR/test_oligos.txt" | tr -d ' ')
    if [ "$ACTUAL_COUNT" -eq 6 ]; then
        # Compare sequences
        awk -F'\t' '{print $(NF-1)}' "$OUTPUT_DIR/test_oligos.txt" | sort > /tmp/val_krt19_actual.txt
        awk -F'\t' '{print $(NF-1)}' test_cases/KRT19_withUTRs/KRT19_withUTRs_oligos.txt | sort > /tmp/val_krt19_expected.txt

        if diff -q /tmp/val_krt19_actual.txt /tmp/val_krt19_expected.txt > /dev/null 2>&1; then
            echo -e "${GREEN}✓ PASS: 6/6 probes match (100%)${NC}"
            echo "Result: PASS - 6/6 probes match" >> "$LOG_FILE"
            ((PASS_COUNT++))
        else
            MATCHING=$(comm -12 /tmp/val_krt19_actual.txt /tmp/val_krt19_expected.txt | wc -l | tr -d ' ')
            echo -e "${RED}✗ FAIL: Only $MATCHING/6 probes match${NC}"
            echo "Result: FAIL - $MATCHING/6 probes match" >> "$LOG_FILE"
            ((FAIL_COUNT++))
        fi
    else
        echo -e "${RED}✗ FAIL: Expected 6 probes, got $ACTUAL_COUNT${NC}"
        echo "Result: FAIL - Expected 6 probes, got $ACTUAL_COUNT" >> "$LOG_FILE"
        ((FAIL_COUNT++))
    fi
else
    echo -e "${RED}✗ FAIL: probedesign command failed${NC}"
    echo "Result: FAIL - Command error" >> "$LOG_FILE"
    ((FAIL_COUNT++))
fi
echo "" >> "$LOG_FILE"
echo ""

# Check 5: EIF1 HCR test
echo -e "${BLUE}[5/5] Running EIF1 HCR test (≥15/20 probes expected)...${NC}"
echo "[5/5] EIF1 HCR test" >> "$LOG_FILE"

OUTPUT_DIR="/tmp/validation_eif1"
rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

if probedesign design test_cases/EIF1_CDS_HCR/EIF1_Exons.fasta \
  --probes 20 \
  --oligo-length 52 \
  --target-gibbs -60 \
  --allowable-gibbs -80,-40 \
  --pseudogene-mask \
  --genome-mask \
  --index-dir bowtie_indexes \
  -o "$OUTPUT_DIR/test" > /dev/null 2>> "$LOG_FILE"; then

    ACTUAL_COUNT=$(wc -l < "$OUTPUT_DIR/test_oligos.txt" | tr -d ' ')

    # Compare sequences
    awk -F'\t' '{print $(NF-1)}' "$OUTPUT_DIR/test_oligos.txt" | sort > /tmp/val_eif1_actual.txt
    awk -F'\t' '{print $(NF-1)}' test_cases/EIF1_CDS_HCR/EIF1_CDS_HCR_oligos.txt | sort > /tmp/val_eif1_expected.txt

    MATCHING=$(comm -12 /tmp/val_eif1_actual.txt /tmp/val_eif1_expected.txt | wc -l | tr -d ' ')
    PERCENT=$((MATCHING * 100 / 20))

    if [ "$MATCHING" -ge 15 ]; then
        echo -e "${GREEN}✓ PASS: $MATCHING/20 probes match ($PERCENT%, ≥75% required)${NC}"
        echo "Result: PASS - $MATCHING/20 probes match ($PERCENT%)" >> "$LOG_FILE"
        ((PASS_COUNT++))
    else
        echo -e "${RED}✗ FAIL: Only $MATCHING/20 probes match ($PERCENT%, <75%)${NC}"
        echo "Result: FAIL - $MATCHING/20 probes match ($PERCENT%)" >> "$LOG_FILE"
        ((FAIL_COUNT++))
    fi
else
    echo -e "${RED}✗ FAIL: probedesign command failed${NC}"
    echo "Result: FAIL - Command error" >> "$LOG_FILE"
    ((FAIL_COUNT++))
fi
echo "" >> "$LOG_FILE"
echo ""

# Summary
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Validation Summary${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo "Total checks: $TOTAL_COUNT"
echo -e "Passed: ${GREEN}$PASS_COUNT${NC}"
echo -e "Failed: ${RED}$FAIL_COUNT${NC}"
echo ""

echo "" >> "$LOG_FILE"
echo "=== Summary ===" >> "$LOG_FILE"
echo "Total: $TOTAL_COUNT" >> "$LOG_FILE"
echo "Passed: $PASS_COUNT" >> "$LOG_FILE"
echo "Failed: $FAIL_COUNT" >> "$LOG_FILE"

if [ "$FAIL_COUNT" -eq 0 ]; then
    echo -e "${GREEN}All validation checks passed! ✓${NC}"
    echo "Status: ALL PASSED" >> "$LOG_FILE"
    exit 0
else
    echo -e "${RED}Some validation checks failed.${NC}"
    echo "Status: SOME FAILED" >> "$LOG_FILE"
    exit 1
fi
