#!/bin/bash
# Download Pre-built Genome Indices
# Purpose: Download pre-built Bowtie 2 indices from AWS (compatible with Bowtie 1)
# Source: https://benlangmead.github.io/aws-indexes/bowtie
# Note: Bowtie 1 can read Bowtie 2 .bt2 index files directly
# Requirements: curl, unzip
# Output: bowtie_indexes/<genome_prefix>.*.bt2 (6 files per genome)

set -euo pipefail

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Paths
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
INDEX_DIR="$REPO_ROOT/bowtie_indexes"
mkdir -p "$INDEX_DIR"

# Log file
LOG_FILE="$REPO_ROOT/jyl/setup/logs/02_genome_download.log"
echo "=== Genome Index Download Log ===" > "$LOG_FILE"
echo "Date: $(date)" >> "$LOG_FILE"
echo "Source: Pre-built Bowtie 1 indices from AWS" >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Downloading Genome Indices${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo "Repository: $REPO_ROOT"
echo "Index directory: $INDEX_DIR"
echo ""

# Function to download and rename index
download_genome() {
    local species="$1"
    local url="$2"
    local zip_name="$3"
    local source_prefix="$4"
    local target_prefix="$5"

    echo -e "${BLUE}=== $species ===${NC}"
    echo "=== $species ===" >> "$LOG_FILE"

    cd "$INDEX_DIR"

    # Check if already exists
    if ls "${target_prefix}".*.bt2 1> /dev/null 2>&1; then
        local count=$(ls "${target_prefix}".*.bt2 | wc -l | tr -d ' ')
        if [ "$count" -eq 6 ]; then
            echo -e "${YELLOW}Already exists (6 files), skipping${NC}"
            echo "Status: Already exists" >> "$LOG_FILE"
            echo "" >> "$LOG_FILE"
            return 0
        fi
    fi

    # Download
    echo "Downloading $zip_name..."
    echo "URL: $url" >> "$LOG_FILE"
    local start=$(date +%s)

    if curl -# -L -o "$zip_name" "$url" >> "$LOG_FILE" 2>&1; then
        local elapsed=$(($(date +%s) - start))
        local size=$(du -h "$zip_name" | cut -f1)
        echo -e "${GREEN}Downloaded: $size in ${elapsed}s${NC}"
        echo "Downloaded: $size in ${elapsed}s" >> "$LOG_FILE"
    else
        echo -e "${RED}ERROR: Download failed${NC}"
        echo "ERROR: Download failed" >> "$LOG_FILE"
        return 1
    fi

    # Extract
    echo "Extracting..."
    start=$(date +%s)
    if unzip -q "$zip_name" >> "$LOG_FILE" 2>&1; then
        elapsed=$(($(date +%s) - start))
        echo -e "${GREEN}Extracted in ${elapsed}s${NC}"
        echo "Extracted in ${elapsed}s" >> "$LOG_FILE"
    else
        echo -e "${RED}ERROR: Extraction failed${NC}"
        echo "ERROR: Extraction failed" >> "$LOG_FILE"
        return 1
    fi

    # Move files from subdirectories if needed
    if [ -d "$source_prefix" ]; then
        echo "Moving files from subdirectory..."
        mv "$source_prefix"/*.bt2 . 2>/dev/null || true
        rmdir "$source_prefix" 2>/dev/null || rm -rf "$source_prefix"
    fi

    # Rename if needed
    if [ "$source_prefix" != "$target_prefix" ]; then
        echo "Renaming $source_prefix → $target_prefix..."
        for file in "${source_prefix}".*.bt2; do
            if [ -f "$file" ]; then
                new_name="${file/$source_prefix/$target_prefix}"
                mv "$file" "$new_name"
            fi
        done
        echo "Renamed to match pipeline expectation" >> "$LOG_FILE"
    fi

    # Verify (.bt2 files work with Bowtie 1)
    local count=$(ls "${target_prefix}".*.bt2 2>/dev/null | wc -l | tr -d ' ')
    if [ "$count" -eq 6 ]; then
        echo -e "${GREEN}Verified: 6 .bt2 files${NC}"
        du -ch "${target_prefix}".*.bt2 | tail -1 | awk '{print "Total size: " $1}'
        echo "Verified: 6 files" >> "$LOG_FILE"
        du -ch "${target_prefix}".*.bt2 >> "$LOG_FILE" 2>&1
    else
        echo -e "${RED}WARNING: Expected 6 files, found $count${NC}"
        echo "WARNING: Found $count files" >> "$LOG_FILE"
    fi

    # Cleanup
    rm -f "$zip_name"
    echo "" >> "$LOG_FILE"
    cd "$REPO_ROOT"
}

# Download indices
echo -e "${BLUE}Downloading pre-built Bowtie 1 indices from AWS...${NC}"
echo ""

# Human GRCh38
download_genome "Human (GRCh38)" \
    "https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip" \
    "GRCh38_noalt_as.zip" \
    "GRCh38_noalt_as" \
    "GCA_000001405.15_GRCh38_no_alt_analysis_set"

# Mouse mm10
download_genome "Mouse (mm10)" \
    "https://genome-idx.s3.amazonaws.com/bt/mm10.zip" \
    "mm10.zip" \
    "mm10" \
    "mm10"

# Drosophila BDGP6
download_genome "Drosophila (BDGP6)" \
    "https://genome-idx.s3.amazonaws.com/bt/BDGP6.zip" \
    "BDGP6.zip" \
    "BDGP6" \
    "drosophila"

# Summary
echo ""
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Download Complete${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo "All genome indices downloaded successfully!"
echo "Location: $INDEX_DIR"
echo ""
echo "Index files:"
ls -1 "$INDEX_DIR"/*.bt2 2>/dev/null | grep -E "(GCA_|mm10|drosophila)" | head -20 | sed 's/^/  /'
echo ""
echo "Total disk usage:"
du -ch "$INDEX_DIR"/*.bt2 2>/dev/null | tail -1 | awk '{print "  " $1}'
echo ""
echo "Log: $LOG_FILE"

# Summary to log
echo "" >> "$LOG_FILE"
echo "=== Summary ===" >> "$LOG_FILE"
echo "All genome indices downloaded" >> "$LOG_FILE"
ls -1 "$INDEX_DIR"/*.bt2 2>/dev/null | grep -E "(GCA_|mm10|drosophila)" >> "$LOG_FILE"
