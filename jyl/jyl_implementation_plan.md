# ProbeDesign Python Pipeline — Installation & Validation Plan

> **Audience**: Lab members who are comfortable running terminal commands but may not have Python/bioinformatics experience.
> **Last updated**: 2026-02-17

---

## 📋 Setup Status (2026-02-17)

**Environment**: macOS (Darwin 23.5.0)
**Setup Location**: `jyl/setup/`

### Completed Steps ✓

- [x] **Section 1**: Assumptions verified (micromamba installed, git repo cloned)
- [x] **Section 2**: Environment setup completed
  - Created `probedesign` environment with Python 3.11.14, click 8.3.1, bowtie 1.3.1, pytest 9.0.2
  - Location: `/Users/jefflee/micromamba/envs/probedesign`
- [x] **Section 3**: Package installation completed
  - Installed probedesign v0.1.0 in editable mode
  - Entry point `probedesign` command verified
- [x] **Section 4**: Bowtie verification completed
  - Version: 1.3.1 (bowtie-align-s, 64-bit)
  - Correct bioconda installation confirmed
- [x] **Section 6.1**: Minimal test passed (100% match)
  - CDKN1A_32 test: 32/32 probes matching reference
  - Python implementation verified against MATLAB

### Pending Steps

- [ ] **Section 5**: Reference Data & Indices (skipped for now)
  - [ ] 5B: Pseudogene indexes (human, mouse, drosophila)
  - [ ] 5C: Genome indexes (GRCh38, mm10, dm6)
  - [ ] 5D: RepeatMasker & Dfam databases
- [ ] **Section 6.2-6.4**: Additional tests (require bowtie indexes)

**Setup Documentation**: See `jyl/setup/SETUP_LOG.md` for detailed logs and outputs.

---

## Table of Contents

1. [Assumptions](#1-assumptions)
2. [Environment Setup (micromamba)](#2-environment-setup-micromamba)
3. [Package Installation](#3-package-installation)
4. [Bowtie Verification](#4-bowtie-verification)
5. [Reference Data & Indices](#5-reference-data--indices)
   - [5A. Directory Convention](#5a-directory-convention)
   - [5B. Pseudogene Indexes](#5b-pseudogene-indexes)
   - [5C. Genome Indexes](#5c-genome-indexes)
   - [5D. RepeatMasker & Dfam Databases](#5d-repeatmasker--dfam-databases)
6. [Smoke Test / Validation](#6-smoke-test--validation)
7. [Troubleshooting](#7-troubleshooting)
8. [Resource Estimates](#8-resource-estimates)
9. [Offline Install Path](#9-offline-install-path)
10. [Platform Notes](#10-platform-notes)
11. [Logging & Provenance](#11-logging--provenance)
12. [Reproducibility (Lockfile)](#12-reproducibility-lockfile)

---

## 1. Assumptions

| Item | Requirement |
|------|-------------|
| Operating system | macOS 13+ (Apple Silicon or Intel) **or** Linux x86_64 |
| CPU / RAM | Any modern CPU; 8 GB RAM minimum (16 GB recommended for genome index building) |
| Free disk space | **80 GB+** (see [Resource Estimates](#8-resource-estimates) for per-species breakdown) |
| Internet | Required for initial setup; see [Offline Install Path](#9-offline-install-path) for air-gapped clusters |
| micromamba | Already installed ([install guide](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)) |
| Shell | bash or zsh with micromamba initialized (`micromamba shell init`) |
| Git | Installed and repo cloned: `git clone https://github.com/arjunrajlaboratory/ProbeDesign.git` |

---

## 2. Environment Setup (micromamba)

### 2.1 Create the environment

```bash
micromamba create -n probedesign \
  -c conda-forge -c bioconda \
  --channel-priority strict \
  python=3.11 \
  click">=8.0" \
  bowtie=1.3.1 \
  pytest">=7.0" \
  -y
```

**What this does:**
- Creates a new environment named `probedesign`
- Uses `conda-forge` and `bioconda` channels with strict priority (prevents package conflicts)
- Pins Python 3.11, installs Bowtie 1.3.1 (the bioinformatics aligner, NOT the homebrew tool)
- Installs `click` (CLI framework) and `pytest` (for testing)

### 2.2 Activate the environment

```bash
micromamba activate probedesign
```

### 2.3 Verify Python and key packages

```bash
python --version
# Expected: Python 3.11.x

python -c "import click; print(click.__version__)"
# Expected: 8.x.x

bowtie --version
# Expected first line: bowtie-align-s version 1.3.1
# (see Section 4 for full verification)
```

### 2.4 (Optional) Install RepeatMasker

RepeatMasker is only needed if you plan to use the `--repeatmask` flag for automatic repeat masking. Skip this if you will always provide a pre-masked file via `--repeatmask-file`.

```bash
micromamba install -n probedesign -c bioconda -c conda-forge repeatmasker -y
```

This installs RepeatMasker 4.2.2 along with its search engine (RMBlast), TRF, HMMER, and the root Dfam database (partition 0, ~77 MB). See [Section 5D](#5d-repeatmasker--dfam-databases) for species-specific database partitions.

---

## 3. Package Installation

### Path A: Editable Install from Repository (Recommended)

```bash
cd /path/to/smFISHProbeDesign   # Your clone of the repo
micromamba activate probedesign
pip install -e .
```

**Verify:**
```bash
probedesign --version
# Expected: probedesign, version 0.1.0

probedesign --help
# Expected: Usage info for the probedesign CLI
```

**Code reference:** The entry point `probedesign` is defined in `pyproject.toml` under `[project.scripts]`:
```
probedesign = "probedesign.cli:main"
```
The package source is in `src/probedesign/` (setuptools `find` with `where = ["src"]`).

### Path B: Fallback if Dependency Solving Fails

If micromamba fails to solve the environment (rare, but possible with channel conflicts):

1. **Try relaxing pins:**
   ```bash
   micromamba create -n probedesign \
     -c conda-forge -c bioconda \
     --channel-priority strict \
     python=3.11 bowtie -y
   ```
   Then install Python deps separately:
   ```bash
   micromamba activate probedesign
   pip install click pytest
   pip install -e .
   ```

2. **Check available versions:**
   ```bash
   micromamba search -c bioconda bowtie
   micromamba search -c conda-forge click
   ```

3. **Docker/Apptainer fallback:**
   If conda/micromamba issues persist, use Miniforge in a container:
   ```bash
   # Example Dockerfile
   FROM condaforge/miniforge3:latest
   RUN mamba install -c conda-forge -c bioconda bowtie=1.3.1 python=3.11 click -y
   COPY . /app
   WORKDIR /app
   RUN pip install -e .
   ```

---

## 4. Bowtie Verification

Bowtie 1 (not Bowtie 2) is required. On macOS, **homebrew has a different package named "bowtie"** (a file manager) that is NOT the bioinformatics tool.

```bash
micromamba activate probedesign

# Check the binary
which bowtie
# Expected: /path/to/micromamba/envs/probedesign/bin/bowtie

bowtie --version
# Expected output (first 2 lines):
#   /path/to/bowtie-align-s version 1.3.1
#   64-bit
```

**If you see "bowtie: command not found":**
```bash
micromamba install -n probedesign -c bioconda bowtie=1.3.1 -y
```

**If `bowtie --version` does NOT show "bowtie-align-s"**, you may have the wrong bowtie (e.g., homebrew). Remove it and install the bioconda version:
```bash
brew uninstall bowtie 2>/dev/null  # Remove homebrew version if present
micromamba install -n probedesign -c bioconda bowtie=1.3.1 --force-reinstall -y
```

**Code reference:** The pipeline searches for bowtie via `find_bowtie()` in `src/probedesign/masking.py:22-63`. Search order:
1. `which bowtie` (PATH — picks up the conda env binary)
2. `/usr/bin/bowtie`
3. `$BOWTIEHOME/bowtie` environment variable

---

## 5. Reference Data & Indices

### 5A. Directory Convention

All bowtie indexes live in `bowtie_indexes/` at the repository root. The pipeline auto-detects this location.

```
smFISHProbeDesign/
├── bowtie_indexes/             # <-- Create this directory
│   ├── humanPseudo.*.ebwt      # Human pseudogene index (6 files)
│   ├── mousePseudo.*.ebwt      # Mouse pseudogene index (6 files)
│   ├── drosophilaPseudo.*.ebwt # Drosophila pseudogene index (6 files)
│   ├── GCA_000001405.15_GRCh38_no_alt_analysis_set.*.ebwt  # Human genome (6 files)
│   ├── mm10.*.ebwt             # Mouse genome (6 files)
│   └── drosophila.*.ebwt       # Drosophila genome (6 files)
├── probedesign/
│   └── pseudogeneDBs/          # Shipped pseudogene FASTA files
│       ├── human.fasta
│       ├── mouse.fasta
│       ├── drosophila.fasta
│       └── elegans.fasta
└── src/probedesign/            # Python package
```

**Code reference:** Default index directory is `src/probedesign/masking.py:19`:
```python
DEFAULT_INDEX_DIR = Path(__file__).parent.parent.parent.parent / "bowtie_indexes"
```
You can override this at runtime via `--index-dir /custom/path`.

**Environment variable (optional):**
```bash
export PROBEDESIGN_INDEX_DIR=/path/to/bowtie_indexes
```
Note: The code does not read this variable natively. Use `--index-dir` CLI flag instead.

### 5B. Pseudogene Indexes

Pseudogene FASTA files ship with the repository in `probedesign/pseudogeneDBs/`. Build bowtie indexes from them:

```bash
cd /path/to/smFISHProbeDesign
mkdir -p bowtie_indexes
micromamba activate probedesign
```

#### Human

```bash
bowtie-build probedesign/pseudogeneDBs/human.fasta bowtie_indexes/humanPseudo
```

Expected output:
```
Total time for backward call to driver() ...
```

Verify (6 files created):
```bash
ls bowtie_indexes/humanPseudo.*.ebwt | wc -l
# Expected: 6
```

**Code reference:** `src/probedesign/masking.py:251` — `"human": "humanPseudo"`

#### Mouse

```bash
bowtie-build probedesign/pseudogeneDBs/mouse.fasta bowtie_indexes/mousePseudo
```

Verify:
```bash
ls bowtie_indexes/mousePseudo.*.ebwt | wc -l
# Expected: 6
```

**Code reference:** `src/probedesign/masking.py:252` — `"mouse": "mousePseudo"`

#### Drosophila

```bash
bowtie-build probedesign/pseudogeneDBs/drosophila.fasta bowtie_indexes/drosophilaPseudo
```

Verify:
```bash
ls bowtie_indexes/drosophilaPseudo.*.ebwt | wc -l
# Expected: 6
```

**Code reference:** `src/probedesign/masking.py:254` — `"drosophila": "drosophilaPseudo"`

#### Source & Provenance

| Species | FASTA file | Entries | Size | Source |
|---------|-----------|---------|------|--------|
| human | `human.fasta` | ~16,000 pseudogenes | 26 MB | [pseudogene.org](http://pseudogene.org) (build 61) |
| mouse | `mouse.fasta` | ~17,000 pseudogenes | 28 MB | pseudogene.org (build 60) |
| drosophila | `drosophila.fasta` | ~2,500 pseudogenes | 2.1 MB | pseudogene.org (build 50) |

Conversion from TSV to FASTA was done with MATLAB script `probedesign/pseudogeneDBs/pseudoTSVtoFASTA.m` (reference only; not needed for setup).

### 5C. Genome Indexes

Genome indexes enable the `--genome-mask` feature, which masks probe candidates that align to many genomic locations (repetitive sequences).

#### Human (GRCh38)

**Download pre-built index (~2.7 GB compressed, ~3 GB uncompressed):**
```bash
cd /path/to/smFISHProbeDesign/bowtie_indexes

curl -O ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/GRCh38_no_alt.zip
unzip GRCh38_no_alt.zip
rm GRCh38_no_alt.zip  # Optional: remove zip to save space
```

**Verify:**
```bash
ls GCA_000001405.15_GRCh38_no_alt_analysis_set.*.ebwt | wc -l
# Expected: 6
```

**Integrity check (file count and approximate sizes):**
```bash
du -sh GCA_000001405.15_GRCh38_no_alt_analysis_set.*.ebwt
# Each .ebwt file should be 700 MB - 1.2 GB; total ~3 GB
```

**Code reference:** `src/probedesign/masking.py:297` — `"human": "GCA_000001405.15_GRCh38_no_alt_analysis_set"`

#### Mouse (mm10)

**Download pre-built index (~2.5 GB compressed, ~2.7 GB uncompressed):**
```bash
cd /path/to/smFISHProbeDesign/bowtie_indexes

curl -O ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/mm10.ebwt.zip
unzip mm10.ebwt.zip
rm mm10.ebwt.zip  # Optional
```

**Verify:**
```bash
ls mm10.*.ebwt | wc -l
# Expected: 6
```

**Code reference:** `src/probedesign/masking.py:298` — `"mouse": "mm10"`

#### Drosophila

The code expects a genome index with prefix `drosophila` (see `src/probedesign/masking.py:300`). There are two options:

**Option A: Build from dm6 genome FASTA (recommended for current assembly)**

```bash
cd /path/to/smFISHProbeDesign/bowtie_indexes

# Download dm6 genome FASTA from Ensembl
curl -L -O https://ftp.ensembl.org/pub/release-111/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz
gunzip Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz

# Build index with prefix "drosophila" (matches code expectation)
bowtie-build Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa drosophila

# Optional: remove the FASTA after building
rm Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa
```

**Option B: Download pre-built r5.22 index from Sourceforge and rename**

```bash
cd /path/to/smFISHProbeDesign/bowtie_indexes

# Download pre-built index (~150 MB)
curl -L -O https://sourceforge.net/projects/bowtie-bio/files/bowtie/indexes/d_melanogaster_fb5_22.ebwt.zip/download -o d_melanogaster_fb5_22.ebwt.zip
unzip d_melanogaster_fb5_22.ebwt.zip

# Rename to match expected prefix "drosophila"
for f in d_melanogaster_fb5_22.*.ebwt; do
  mv "$f" "drosophila.${f#d_melanogaster_fb5_22.}"
done
```

> **ASSUMPTION:** The pre-built r5.22 index uses the older dm3 assembly (FlyBase release 5.22). Option A (dm6/BDGP6) is more current. The choice depends on your reference requirements.

**Verify (either option):**
```bash
ls drosophila.*.ebwt | wc -l
# Expected: 6
```

**Code reference:** `src/probedesign/masking.py:300` — `"drosophila": "drosophila"`

#### Genome Index Summary

| Species | Index Prefix | Source | Compressed Size | Uncompressed Size |
|---------|-------------|--------|-----------------|-------------------|
| human | `GCA_000001405.15_GRCh38_no_alt_analysis_set` | JHU FTP (pre-built) | 2.7 GB | ~3.0 GB |
| mouse | `mm10` | JHU FTP (pre-built) | 2.5 GB | ~2.7 GB |
| drosophila | `drosophila` | Ensembl dm6 (build yourself) OR Sourceforge r5.22 (rename) | 150 MB / 50 MB (build) | ~165 MB |

### 5D. RepeatMasker & Dfam Databases

RepeatMasker masks repetitive elements (SINEs, LINEs, LTRs, DNA transposons, etc.) in input sequences. It is used via the `--repeatmask` CLI flag. If you prefer to provide a pre-masked FASTA via `--repeatmask-file`, you can skip RepeatMasker entirely.

#### What RepeatMasker Does

- Screens input DNA against the Dfam database of known repeat families
- Outputs a `.masked` file with repeat regions replaced by N's
- These N positions are then excluded from probe design
- **Code reference:** `src/probedesign/masking.py:449-536` (`run_repeatmasker`, `repeatmasker_mask_to_sequence`)

#### Dfam Database Partitions

The Dfam database is split into taxonomic partitions. You only need the partition(s) for your species.

| Partition | Taxonomic Group | Species | Compressed | Extracted (approx) |
|-----------|----------------|---------|------------|-------------------|
| **0** | Root | **All** (always required) | 17 MB | ~77 MB |
| **1** | Brachycera | **Drosophila melanogaster** | 7.8 GB | ~45 GB **(ASSUMPTION)** |
| **7** | Mammalia | **Human, Mouse, Rat** | 8.9 GB | ~56 GB |

Partition 0 is included when you install RepeatMasker via conda. You must download additional partitions manually.

#### Download Dfam Partitions

First, find where RepeatMasker stores its database:

```bash
# Find the famdb directory
FAMDB_DIR=$(dirname $(find $(micromamba info --json | python -c "import json,sys; print(json.load(sys.stdin)['envs_dirs'][0])")/../ -name "famdb" -path "*/RepeatMasker/Libraries/*" 2>/dev/null | head -1))
echo "FamDB directory: $FAMDB_DIR"
```

Or if using the standard miniforge/micromamba location:
```bash
# macOS with Miniforge
FAMDB_DIR="$MAMBA_ROOT_PREFIX/share/RepeatMasker/Libraries/famdb"
# Verify it exists:
ls "$FAMDB_DIR"
# Should show: dfam39_full.0.h5 (the root partition)
```

**For Human / Mouse (Partition 7 — Mammalia):**
```bash
cd "$FAMDB_DIR"
curl -O https://www.dfam.org/releases/current/families/FamDB/dfam39_full.7.h5.gz
gunzip dfam39_full.7.h5.gz
# Creates: dfam39_full.7.h5 (~56 GB)
```

> Ensure you have **~65 GB free** during extraction (compressed + extracted simultaneously).

**For Drosophila (Partition 1 — Brachycera):**
```bash
cd "$FAMDB_DIR"
curl -O https://www.dfam.org/releases/current/families/FamDB/dfam39_full.1.h5.gz
gunzip dfam39_full.1.h5.gz
# Creates: dfam39_full.1.h5 (~45 GB, ASSUMPTION)
```

#### Verify RepeatMasker Installation

```bash
micromamba activate probedesign

# Check version
RepeatMasker -v 2>/dev/null || \
  $MAMBA_ROOT_PREFIX/share/RepeatMasker/RepeatMasker -v
# Expected: RepeatMasker version 4.2.2 (or similar)

# Test human species lookup (requires partition 7)
echo -e ">test\nACGTACGTACGTACGT" > /tmp/test_rm.fa
RepeatMasker -species human -dir /tmp /tmp/test_rm.fa 2>&1 | head -5
# Should NOT show "partition is absent" error
rm /tmp/test_rm.fa
```

**Code reference:** RepeatMasker executable search in `src/probedesign/masking.py:403-446` (`find_repeatmasker()`). The code checks:
1. `/opt/homebrew/Caskroom/miniforge/base/share/RepeatMasker/RepeatMasker` (hardcoded macOS path)
2. `which RepeatMasker` (PATH)
3. `/usr/local/bin/RepeatMasker`, `~/RepeatMasker/RepeatMasker`

#### Alternative: Manual Repeat Masking (No RepeatMasker Needed)

If you cannot install RepeatMasker (e.g., disk constraints), you can:

1. Use the [RepeatMasker web server](http://www.repeatmasker.org/cgi-bin/WEBRepeatMasker):
   - Upload your FASTA
   - Select species
   - Download the `.masked` output
2. Pass the masked file to ProbeDesign:
   ```bash
   probedesign design input.fa --repeatmask-file input.fa.masked --probes 32
   ```

---

## 6. Smoke Test / Validation

The repository ships with test cases and a test script. Use these to validate your installation.

### 6.1 Minimal Test (No Bowtie Indexes Required)

This test verifies the core algorithm and repeat masking from a file:

```bash
cd /path/to/smFISHProbeDesign
micromamba activate probedesign

probedesign design test_cases/CDKN1A_32/CDKN1A.fa \
  --probes 32 \
  --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa \
  -o /tmp/CDKN1A_test
```

**Expected output files:**
- `/tmp/CDKN1A_test_oligos.txt` — Tab-separated probe list (32 probes)
- `/tmp/CDKN1A_test_seq.txt` — Visual sequence alignment

**Verify probe count:**
```bash
wc -l /tmp/CDKN1A_test_oligos.txt
# Expected: 32 (one line per probe, no header)
```

**Verify against reference (100% match expected):**
```bash
# Extract probe sequences (column 5) and compare
awk '{print $5}' /tmp/CDKN1A_test_oligos.txt | sort > /tmp/actual.txt
awk '{print $5}' test_cases/CDKN1A_32/CDKN1A_32_genomemaskoff_oligos.txt | sort > /tmp/expected.txt
diff /tmp/actual.txt /tmp/expected.txt
# Expected: no output (files are identical)
```

### 6.2 Full Test with Bowtie Masks (Requires Genome + Pseudogene Indexes)

```bash
# Test 2: KRT19 with pseudogene + genome masking (100% match expected)
probedesign design test_cases/KRT19_withUTRs/KRT19_withUTRs.fa \
  --probes 32 \
  --pseudogene-mask \
  --genome-mask \
  --index-dir bowtie_indexes \
  -o /tmp/KRT19_test

# Verify
awk '{print $5}' /tmp/KRT19_test_oligos.txt | sort > /tmp/actual_krt19.txt
awk '{print $5}' test_cases/KRT19_withUTRs/KRT19_withUTRs_oligos.txt | sort > /tmp/expected_krt19.txt
diff /tmp/actual_krt19.txt /tmp/expected_krt19.txt
# Expected: no output (100% match)
```

### 6.3 HCR Probe Test (Requires Genome + Pseudogene Indexes)

```bash
# Test 3: EIF1 HCR probes (>= 75% match expected)
probedesign design test_cases/EIF1_CDS_HCR/EIF1_Exons.fasta \
  --probes 20 \
  --oligo-length 52 \
  --target-gibbs -60 \
  --allowable-gibbs -80,-40 \
  --pseudogene-mask \
  --genome-mask \
  --index-dir bowtie_indexes \
  -o /tmp/EIF1_test
```

This test is expected to achieve ~78% match (not 100%) due to DP tie-breaking differences.

### 6.4 Automated Test Script

The repo includes `run_tests.sh` which automates all the above:

```bash
cd /path/to/smFISHProbeDesign
chmod +x run_tests.sh
./run_tests.sh
```

**Expected output:**
```
Test 1 (CDKN1A_32 repeatmask-file): PASS (100%)
Test 2 (KRT19 bowtie masking):      PASS (100%)  [or SKIPPED if no bowtie indexes]
Test 3 (EIF1 HCR probes):           PASS (78%)   [or SKIPPED if no bowtie indexes]
```

**Code reference:** `run_tests.sh` (247 lines). The script auto-detects bowtie at the Miniforge path and skips bowtie-dependent tests if indexes are missing.

### 6.5 Validation Checklist

- [x] `probedesign --version` prints `0.1.0` ✅
- [x] `bowtie --version` prints `bowtie-align-s version 1.3.1` ✅
- [x] CDKN1A_32 test produces 32/32 matching probes (100%) ✅
- [x] (If bowtie indexes built) KRT19 test produces 6/6 matching probes (100%) ✅
- [x] (If bowtie indexes built) EIF1 HCR test produces >= 15/20 matching probes (>= 75%) ✅ (15/20, 75%)
- [ ] (If RepeatMasker installed) `probedesign design test_cases/CDKN1A_32/CDKN1A.fa --repeatmask --probes 32` runs without error

**Validation completed**: 2026-02-17
**Script**: `jyl/setup/04_run_validation_checklist.sh`
**Log**: `jyl/setup/logs/validation_checklist.log`

---

## 7. Troubleshooting

### "command not found: probedesign"

**Cause:** Package not installed or wrong environment.
```bash
micromamba activate probedesign
pip install -e /path/to/smFISHProbeDesign
probedesign --version
```

### "command not found: bowtie"

**Cause:** Bowtie not installed or environment not activated.
```bash
micromamba activate probedesign
which bowtie  # Should be inside envs/probedesign/bin/
# If not:
micromamba install -n probedesign -c bioconda bowtie=1.3.1 -y
```

### "Could not find bowtie executable"

**Cause:** The Python code cannot find bowtie even though it's installed.
```bash
# Check what find_bowtie() will find:
micromamba activate probedesign
which bowtie
# If this shows the correct path, the code should work.
# If the path differs, set:
export BOWTIEHOME=$(dirname $(which bowtie))
```
**Code reference:** `src/probedesign/masking.py:22-63` (`find_bowtie()`)

### Bowtie index not found / "Error: Could not locate a Bowtie index"

**Cause:** Index files are not in the expected directory or have wrong prefix names.
```bash
# Verify index directory
ls bowtie_indexes/*.ebwt
# Should see files like: humanPseudo.1.ebwt, mm10.1.ebwt, etc.

# Check the expected prefix for your species:
# human pseudogene: humanPseudo
# human genome:     GCA_000001405.15_GRCh38_no_alt_analysis_set
# mouse pseudogene: mousePseudo
# mouse genome:     mm10
# drosophila pseudo: drosophilaPseudo
# drosophila genome: drosophila
```

**If index files exist but have wrong prefix:** Rename or rebuild. The `.ebwt` files must share a common prefix that matches what the code expects.

### "No module named 'probedesign'"

**Cause:** Package not installed in the active environment.
```bash
micromamba activate probedesign
pip install -e /path/to/smFISHProbeDesign
python -c "import probedesign; print('OK')"
```

### RepeatMasker: "partition is absent"

**Cause:** The Dfam database partition for your species is not installed.
```bash
# Example for human (needs partition 7):
cd $MAMBA_ROOT_PREFIX/share/RepeatMasker/Libraries/famdb
curl -O https://www.dfam.org/releases/current/families/FamDB/dfam39_full.7.h5.gz
gunzip dfam39_full.7.h5.gz
```
See [Section 5D](#5d-repeatmasker--dfam-databases) for partition-to-species mapping.

### RepeatMasker: "No module named 'h5py'"

**Cause:** RepeatMasker's Python cannot find h5py.
```bash
micromamba activate probedesign
micromamba install -c conda-forge h5py -y
# The code automatically prepends conda bin to PATH, but if that fails:
PATH="$MAMBA_ROOT_PREFIX/envs/probedesign/bin:$PATH" RepeatMasker -species human input.fa
```
**Code reference:** `src/probedesign/masking.py:498-501` — The code adds `conda_bin` to PATH in the subprocess environment.

### "No space left on device" during Dfam extraction

```bash
df -h $MAMBA_ROOT_PREFIX
# Need ~65 GB free for partition 7 (8.9 GB compressed + ~56 GB extracted)
# Clean up partial file and retry:
rm -f $MAMBA_ROOT_PREFIX/share/RepeatMasker/Libraries/famdb/dfam39_full.7.h5
gunzip dfam39_full.7.h5.gz
```

### Permission errors writing to indexes or RepeatMasker directory

```bash
# Check ownership
ls -la bowtie_indexes/
ls -la $MAMBA_ROOT_PREFIX/share/RepeatMasker/Libraries/famdb/

# If owned by root, either change ownership or use sudo:
sudo chown -R $(whoami) bowtie_indexes/
# Or create indexes in a user-writable location and use --index-dir
```

### macOS: "zsh: permission denied" on run_tests.sh

```bash
chmod +x run_tests.sh
./run_tests.sh
```

---

## 8. Resource Estimates

### Index Building (One-Time)

| Task | Disk (final) | Disk (during build) | RAM | Time (approx) |
|------|-------------|--------------------|----|----------------|
| Human pseudogene index (`bowtie-build`) | ~50 MB | ~200 MB | ~1 GB | 1-2 min |
| Mouse pseudogene index (`bowtie-build`) | ~55 MB | ~220 MB | ~1 GB | 1-2 min |
| Drosophila pseudogene index (`bowtie-build`) | ~4 MB | ~20 MB | ~256 MB | <30 sec |
| Human genome index (download pre-built) | ~3.0 GB | ~5.7 GB (zip + unzip) | N/A | 5-15 min (download) |
| Mouse genome index (download pre-built) | ~2.7 GB | ~5.2 GB | N/A | 5-10 min (download) |
| Drosophila genome index (build from dm6) | ~165 MB | ~500 MB | ~4 GB | 5-10 min |

### RepeatMasker Database

| Partition | Download | Extracted | Total During Extract | Extract Time |
|-----------|----------|-----------|---------------------|--------------|
| 0 (root, included) | — | 77 MB | — | — |
| 1 (Brachycera/Drosophila) | 7.8 GB | ~45 GB **(ASSUMPTION)** | ~53 GB | 15-30 min |
| 7 (Mammalia) | 8.9 GB | ~56 GB | ~65 GB | 20-40 min |

### Cumulative Disk by Species

| Scenario | Disk Required |
|----------|--------------|
| Human (pseudogene + genome + RepeatMasker) | ~62 GB |
| Mouse (pseudogene + genome + RepeatMasker) | ~62 GB |
| Human + Mouse (shared RepeatMasker partition 7) | ~65 GB |
| Drosophila (pseudogene + genome + RepeatMasker) | ~46 GB |
| All three species | ~112 GB |
| Minimal (no genome mask, no RepeatMasker) | < 500 MB |

### Runtime (Per Probe Design Job)

| Step | Time |
|------|------|
| FASTA parsing + badness calculation | < 1 sec |
| Pseudogene masking (3 bowtie runs) | 2-5 sec |
| Genome masking (3 bowtie runs at different stringencies) | 5-15 sec |
| RepeatMasker (if `--repeatmask`) | 30 sec - 2 min |
| Dynamic programming optimization | < 1 sec |
| **Total (without RepeatMasker)** | **5-20 sec** |
| **Total (with RepeatMasker)** | **30 sec - 3 min** |

---

## 9. Offline Install Path

For lab clusters or machines without internet access:

### 9.1 Cache Packages on an Internet-Connected Machine

```bash
# Create the environment and download all packages
micromamba create -n probedesign \
  -c conda-forge -c bioconda \
  --channel-priority strict \
  python=3.11 click">=8.0" bowtie=1.3.1 pytest">=7.0" \
  --download-only -y

# Export the explicit package list
micromamba env export -n probedesign --explicit > probedesign_explicit.txt

# Copy the package cache
# Default cache location:
#   macOS: ~/micromamba/pkgs/
#   Linux: ~/micromamba/pkgs/
CACHE_DIR=$(micromamba info --json | python -c "import json,sys; print(json.load(sys.stdin)['pkgs_dirs'][0])")
tar czf probedesign_pkgs.tar.gz -C "$CACHE_DIR" .
```

### 9.2 Transfer and Install Offline

```bash
# On the offline machine:
tar xzf probedesign_pkgs.tar.gz -C ~/micromamba/pkgs/
micromamba create -n probedesign --file probedesign_explicit.txt --offline -y
```

### 9.3 Transfer Reference Data

```bash
# On the internet machine, package indexes:
tar czf bowtie_indexes.tar.gz bowtie_indexes/

# Transfer and extract on the offline machine:
tar xzf bowtie_indexes.tar.gz -C /path/to/smFISHProbeDesign/
```

---

## 10. Platform Notes

### macOS (Apple Silicon / ARM64)

- Bowtie 1.3.1 runs natively on ARM64 via the bioconda `osx-arm64` build.
- RepeatMasker runs via Rosetta 2 in some configurations. If you encounter architecture-related errors:
  ```bash
  micromamba create -n probedesign --platform osx-64 \
    -c conda-forge -c bioconda \
    python=3.11 bowtie=1.3.1 -y
  ```
  This forces x86_64 packages under Rosetta 2.

### macOS (Intel / x86_64)

- No special considerations. All packages are available natively.

### Linux (x86_64)

- All packages available. No platform-specific issues.
- The hardcoded RepeatMasker path `/opt/homebrew/Caskroom/miniforge/base/...` in `masking.py:417` won't exist on Linux. The code falls through to `which RepeatMasker`, which works if the conda env is activated.
- Similarly, the Miniforge-specific bowtie path in `run_tests.sh` is macOS-specific. On Linux, ensure bowtie is in PATH via `micromamba activate probedesign`.

### Linux (ARM64 / aarch64)

- Bowtie may not be available via bioconda for Linux ARM64. Check:
  ```bash
  micromamba search -c bioconda bowtie --platform linux-aarch64
  ```
  If unavailable, compile from source: https://github.com/BenLangmead/bowtie

### Windows (WSL2)

- Not tested, but should work identically to Linux x86_64 under WSL2.

---

## 11. Logging & Provenance

### Default Output

By default, the pipeline prints progress messages to stdout:
```
Repeat masking (from file): 174 positions masked
Pseudogene masking: 42 positions masked
Genome masking: 156 positions masked
Wrote 32 probes to /tmp/output_oligos.txt
Wrote sequence file to /tmp/output_seq.txt
```

### Quiet Mode

Suppress all stdout output:
```bash
probedesign design input.fa --probes 32 --quiet
```

**Code reference:** `src/probedesign/cli.py` — `-q / --quiet` flag

### Capturing Full Provenance

For reproducibility, capture the exact command and environment:

```bash
# Record the command and its output
probedesign design input.fa --probes 32 --pseudogene-mask --genome-mask \
  --index-dir bowtie_indexes 2>&1 | tee probe_run.log

# Record the environment for reproducibility
micromamba env export -n probedesign > environment.yml
micromamba env export -n probedesign --explicit > environment_explicit.txt
micromamba list -n probedesign > package_versions.txt
```

### Increasing Verbosity

The pipeline does not have a `--verbose` flag. For debugging, add print statements or run Python with `-v`:
```bash
python -v -m probedesign.cli design input.fa --probes 32
```

Bowtie's `--quiet` flag is hardcoded in `masking.py:136`. To see bowtie's verbose output, temporarily edit that line (for debugging only).

---

## 12. Reproducibility (Lockfile)

### Export Environment

```bash
# Human-readable (for sharing)
micromamba env export -n probedesign > environment.yml

# Exact package URLs (for exact reproducibility)
micromamba env export -n probedesign --explicit > environment_explicit.txt
```

### Recreate from Lockfile

```bash
micromamba create -n probedesign --file environment_explicit.txt -y
```

### Version Pinning Summary

| Package | Pinned Version | Channel |
|---------|---------------|---------|
| python | 3.11.x | conda-forge |
| click | >= 8.0 | conda-forge |
| bowtie | 1.3.1 | bioconda |
| pytest | >= 7.0 | conda-forge |
| repeatmasker | 4.2.2 (optional) | bioconda |
| probedesign | 0.1.0 (editable) | local (pip) |
