#!/bin/bash
# Build Bowtie Pseudogene Indices
# Purpose: Build bowtie1 indices from pseudogene FASTA files shipped with the repository
# Requirements: bowtie-build (from probedesign conda environment)
# Output: bowtie_indexes/<species>Pseudo.*.ebwt (6 files per species)

set -euo pipefail

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Get repository root (parent of jyl directory)
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
echo -e "${BLUE}Repository root: $REPO_ROOT${NC}"

# Create bowtie_indexes directory
INDEX_DIR="$REPO_ROOT/bowtie_indexes"
mkdir -p "$INDEX_DIR"
echo -e "${BLUE}Index directory: $INDEX_DIR${NC}"

# Source FASTA directory
FASTA_DIR="$REPO_ROOT/probedesign/pseudogeneDBs"
echo -e "${BLUE}FASTA directory: $FASTA_DIR${NC}"

# Output log file
LOG_FILE="$REPO_ROOT/jyl/pseudogene_index_build.log"
echo "=== Pseudogene Index Build Log ===" > "$LOG_FILE"
echo "Date: $(date)" >> "$LOG_FILE"
echo "Repository: $REPO_ROOT" >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

# Check bowtie-build availability
if ! command -v bowtie-build &> /dev/null; then
    echo -e "${RED}ERROR: bowtie-build not found in PATH${NC}"
    echo "Please activate the probedesign conda environment:"
    echo "  micromamba activate probedesign"
    exit 1
fi

BOWTIE_VERSION=$(bowtie-build --version 2>&1 | head -1)
echo -e "${GREEN}Found: $BOWTIE_VERSION${NC}"
echo "Bowtie: $BOWTIE_VERSION" >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

# Function to build index
build_index() {
    local species="$1"
    local fasta_file="$2"
    local index_prefix="$3"

    echo ""
    echo -e "${BLUE}=== Building $species pseudogene index ===${NC}"
    echo "=== $species ===" >> "$LOG_FILE"
    echo "Input FASTA: $fasta_file" >> "$LOG_FILE"
    echo "Output prefix: $index_prefix" >> "$LOG_FILE"

    if [ ! -f "$fasta_file" ]; then
        echo -e "${RED}ERROR: FASTA file not found: $fasta_file${NC}"
        echo "ERROR: FASTA file not found" >> "$LOG_FILE"
        return 1
    fi

    # Get FASTA stats
    local num_seqs=$(grep -c "^>" "$fasta_file")
    local file_size=$(du -h "$fasta_file" | cut -f1)
    echo "FASTA stats: $num_seqs sequences, $file_size"
    echo "FASTA stats: $num_seqs sequences, $file_size" >> "$LOG_FILE"

    # Build index
    echo "Running bowtie-build..."
    local start_time=$(date +%s)

    if bowtie-build "$fasta_file" "$INDEX_DIR/$index_prefix" >> "$LOG_FILE" 2>&1; then
        local end_time=$(date +%s)
        local elapsed=$((end_time - start_time))
        echo -e "${GREEN}SUCCESS: Built in ${elapsed}s${NC}"
        echo "Build time: ${elapsed}s" >> "$LOG_FILE"

        # Verify output files
        local ebwt_count=$(ls "$INDEX_DIR/$index_prefix".*.ebwt 2>/dev/null | wc -l | tr -d ' ')
        if [ "$ebwt_count" -eq 6 ]; then
            echo -e "${GREEN}Verified: 6 .ebwt files created${NC}"
            echo "Verified: 6 .ebwt files" >> "$LOG_FILE"

            # Show file sizes
            du -h "$INDEX_DIR/$index_prefix".*.ebwt | tail -1 | awk '{print "Total size: " $1}'
            du -ch "$INDEX_DIR/$index_prefix".*.ebwt | tail -1 | awk '{print "Total size: " $1}' >> "$LOG_FILE"
        else
            echo -e "${RED}WARNING: Expected 6 .ebwt files, found $ebwt_count${NC}"
            echo "WARNING: Expected 6 files, found $ebwt_count" >> "$LOG_FILE"
        fi
    else
        echo -e "${RED}ERROR: bowtie-build failed${NC}"
        echo "ERROR: bowtie-build failed" >> "$LOG_FILE"
        return 1
    fi

    echo "" >> "$LOG_FILE"
}

# Build indices for each species
echo ""
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Starting Pseudogene Index Build${NC}"
echo -e "${BLUE}========================================${NC}"

# Human
build_index "human" "$FASTA_DIR/human.fasta" "humanPseudo"

# Mouse
build_index "mouse" "$FASTA_DIR/mouse.fasta" "mousePseudo"

# Drosophila
build_index "drosophila" "$FASTA_DIR/drosophila.fasta" "drosophilaPseudo"

# Summary
echo ""
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Build Complete${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo "All pseudogene indices built successfully!"
echo "Index location: $INDEX_DIR"
echo ""
echo "Index files created:"
ls -1 "$INDEX_DIR"/*Pseudo*.ebwt | sed 's/^/  /'
echo ""
echo "Disk usage:"
du -h "$INDEX_DIR"/*Pseudo*.ebwt | awk '{sum+=$1} END {print "  Total: ~" int(sum) " KB"}'
echo ""
echo "Log file: $LOG_FILE"
echo ""

# Append summary to log
echo "" >> "$LOG_FILE"
echo "=== Summary ===" >> "$LOG_FILE"
echo "All indices built successfully" >> "$LOG_FILE"
echo "Location: $INDEX_DIR" >> "$LOG_FILE"
ls -1 "$INDEX_DIR"/*Pseudo*.ebwt >> "$LOG_FILE"
