#!/bin/bash
# Test Bowtie Indices
# Purpose: Perform minimal queries against pseudogene and genome indices to verify they work
# Requirements: bowtie (from probedesign conda environment), index files in bowtie_indexes/

set -euo pipefail

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get repository root
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
INDEX_DIR="$REPO_ROOT/bowtie_indexes"

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Testing Bowtie Indices${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Check bowtie availability
if ! command -v bowtie &> /dev/null; then
    echo -e "${RED}ERROR: bowtie not found in PATH${NC}"
    echo "Please activate the probedesign conda environment:"
    echo "  micromamba activate probedesign"
    exit 1
fi

BOWTIE_VERSION=$(bowtie --version 2>&1 | head -1)
echo -e "${GREEN}Found: $BOWTIE_VERSION${NC}"
echo ""

# Create test sequence (16bp, simple)
TEST_SEQ="ATCGATCGATCGATCG"
TEST_FASTA="/tmp/test_probe.fa"
echo ">test_probe" > "$TEST_FASTA"
echo "$TEST_SEQ" >> "$TEST_FASTA"

# Function to test an index
test_index() {
    local index_name="$1"
    local index_prefix="$2"

    echo -e "${BLUE}Testing $index_name index...${NC}"

    # Check if index files exist (.ebwt for pseudogenes, .bt2 for genomes)
    local ebwt_count=$(ls "$INDEX_DIR/$index_prefix".*.ebwt 2>/dev/null | wc -l | tr -d ' ')
    local bt2_count=$(ls "$INDEX_DIR/$index_prefix".*.bt2 2>/dev/null | wc -l | tr -d ' ')
    local total_count=$((ebwt_count + bt2_count))

    if [ "$total_count" -ne 6 ]; then
        echo -e "${RED}  SKIP: Index not found or incomplete ($total_count/6 files)${NC}"
        echo ""
        return 1
    fi

    local file_type="Bowtie 1 (.ebwt)"
    if [ "$bt2_count" -eq 6 ]; then
        file_type="Bowtie 2 (.bt2, compatible)"
    fi

    # Run bowtie query (16bp exact match, report up to 10 hits)
    local output_file="/tmp/test_${index_prefix}.txt"
    if bowtie -f -v 0 -k 10 --quiet "$INDEX_DIR/$index_prefix" "$TEST_FASTA" > "$output_file" 2>&1; then
        local hit_count=$(wc -l < "$output_file" | tr -d ' ')
        echo -e "${GREEN}  PASS: Query executed successfully${NC}"
        echo "  Format: $file_type"
        echo "  Hits found: $hit_count (for test sequence: $TEST_SEQ)"
        echo ""
    else
        echo -e "${RED}  FAIL: Bowtie query failed${NC}"
        echo ""
        return 1
    fi
}

# Test pseudogene indices
echo -e "${BLUE}--- Pseudogene Indices ---${NC}"
echo ""
test_index "Human pseudogene" "humanPseudo"
test_index "Mouse pseudogene" "mousePseudo"
test_index "Drosophila pseudogene" "drosophilaPseudo"

# Test genome indices (may not exist yet)
echo -e "${BLUE}--- Genome Indices ---${NC}"
echo ""
test_index "Human genome (GRCh38)" "GCA_000001405.15_GRCh38_no_alt_analysis_set"
test_index "Mouse genome (mm10)" "mm10"
test_index "Drosophila genome (dm6)" "drosophila"

# Summary
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Index Testing Complete${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo "Test FASTA used: $TEST_FASTA"
echo "Test sequence: $TEST_SEQ (16bp)"
echo ""
echo "All available indices have been tested with a minimal query."
echo "Pseudogene indices: Bowtie 1 (.ebwt) format"
echo "Genome indices: Bowtie 2 (.bt2) format (compatible with Bowtie 1)"
echo ""

# Clean up
rm -f "$TEST_FASTA" /tmp/test_*.txt
