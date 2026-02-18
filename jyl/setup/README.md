# ProbeDesign Setup Guide

Complete! All setup steps finished and tested.

**Last Updated**: 2026-02-17
**Location**: `jyl/setup/`

---

## Setup Status: ✅ Complete

| Step | Status | Time | Size |
|------|--------|------|------|
| 1. Environment | ✅ Complete | 5 min | - |
| 2. Pseudogene indices | ✅ Complete | 47 sec | 84 MB |
| 3. Genome indices | ✅ Complete | 12 min | 6.7 GB |
| 4. Testing | ✅ Passed | <1 min | - |

---

## Quick Verification

### Check environment:
```bash
micromamba activate probedesign
probedesign --version  # Should show: 0.1.0
bowtie --version       # Should show: bowtie-align-s version 1.3.1
```

### Test indices:
```bash
cd /Users/jefflee/Documents/Github/smFISHProbeDesign
./jyl/setup/03_test_indices.sh
```

Expected: All 6 indices (3 pseudogene + 3genome) pass.

---

## What Was Built

### Step 1: Environment Setup ✅
- Python 3.11.14
- bowtie 1.3.1
- click, pytest
- probedesign v0.1.0

**Documentation**: [`docs/01_environment_setup.md`](docs/01_environment_setup.md)

---

### Step 2: Pseudogene Indices ✅

**Script**: [`01_build_pseudogene_indices.sh`](01_build_pseudogene_indices.sh)
**Log**: [`logs/01_pseudogene_build.log`](logs/01_pseudogene_build.log)

Built from `probedesign/pseudogeneDBs/*.fasta`:
- **Human**: `humanPseudo` (18,046 sequences, 35 MB, built in 23s)
- **Mouse**: `mousePseudo` (19,086 sequences, 39 MB, built in 22s)
- **Drosophila**: `drosophilaPseudo` (2,204 sequences, 10 MB, built in 2s)

**Format**: Bowtie 1 (`.ebwt`) — 6 files per species, 84 MB total

---

### Step 3: Genome Indices ✅

**Script**: [`02_download_genome_indices.sh`](02_download_genome_indices.sh)
**Log**: [`logs/02_genome_download.log`](logs/02_genome_download.log)

Downloaded pre-built indices from AWS (https://benlangmead.github.io/aws-indexes/bowtie):
- **Human (GRCh38)**: `GCA_000001405.15_GRCh38_no_alt_analysis_set` (3.5 GB, 331s download)
- **Mouse (mm10)**: `mm10` (3.1 GB, 296s download)
- **Drosophila (BDGP6)**: `drosophila` (176 MB, 17s download)

**Format**: Bowtie 2 (`.bt2`) — 6 files per genome, 6.7 GB total
**Note**: Bowtie 1 can read Bowtie 2 indices directly — tested and working!

---

### Step 4: Index Testing ✅

**Script**: [`03_test_indices.sh`](03_test_indices.sh)

All 6 indices tested with minimal bowtie queries:
- ✅ Human pseudogene (5 hits)
- ✅ Mouse pseudogene (7 hits)
- ✅ Drosophila pseudogene (5 hits)
- ✅ Human genome (12 hits)
- ✅ Mouse genome (15 hits)
- ✅ Drosophila genome (15 hits)

---

## Index Locations

```
smFISHProbeDesign/bowtie_indexes/
├── humanPseudo.{1,2,3,4,rev.1,rev.2}.ebwt          # 35 MB
├── mousePseudo.{1,2,3,4,rev.1,rev.2}.ebwt          # 39 MB
├── drosophilaPseudo.{1,2,3,4,rev.1,rev.2}.ebwt     # 10 MB
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.{1,2,3,4,rev.1,rev.2}.bt2    # 3.5 GB
├── mm10.{1,2,3,4,rev.1,rev.2}.bt2                   # 3.1 GB
└── drosophila.{1,2,3,4,rev.1,rev.2}.bt2             # 176 MB

Total: 18 pseudogene files (84 MB) + 18 genome files (6.7 GB) = 6.8 GB
```

---

## Pipeline Usage

### Basic probe design:
```bash
probedesign design input.fa --probes 32
```

### With repeat masking:
```bash
probedesign design input.fa --probes 32 --repeatmask-file masked.fa
```

### With pseudogene masking:
```bash
probedesign design input.fa --probes 32 --pseudogene-mask
```

### With genome masking:
```bash
probedesign design input.fa --probes 32 --genome-mask
```

### Full masking (all features):
```bash
probedesign design input.fa --probes 32 \
  --repeatmask-file masked.fa \
  --pseudogene-mask \
  --genome-mask
```

---

## Validation Tests

Run full pipeline validation:

```bash
cd /Users/jefflee/Documents/Github/smFISHProbeDesign
./run_tests.sh
```

**Expected results**:
- Test 1 (CDKN1A_32): 100% match ✅
- Test 2 (KRT19): 100% match ✅ (requires genome indices)
- Test 3 (EIF1_HCR): 78% match ✅ (requires genome indices)

### Validation Checklist

Run automated validation checklist (Section 6.5 of implementation plan):

```bash
./jyl/setup/04_run_validation_checklist.sh
```

**Results** (2026-02-17):
- ✅ probedesign version 0.1.0
- ✅ bowtie version 1.3.1
- ✅ CDKN1A_32: 32/32 probes (100%)
- ✅ KRT19: 6/6 probes (100%)
- ✅ EIF1 HCR: 15/20 probes (75%)

**All 5 validation checks passed!**

---

## Directory Structure

```
jyl/setup/
├── README.md                           # This file
│
├── 01_build_pseudogene_indices.sh      # Builds pseudogene indices ✅
├── 02_download_genome_indices.sh       # Downloads genome indices ✅
├── 03_test_indices.sh                  # Tests all indices ✅
├── 04_run_validation_checklist.sh      # Run validation checklist ✅
│
├── docs/                               # Documentation
│   ├── 00_index_setup_guide.md         # Detailed index setup guide
│   ├── 01_environment_setup.md         # Environment setup docs
│   └── 01_environment_setup_log.md     # Setup log
│
├── logs/                               # Build logs
│   ├── 01_pseudogene_build.log         # Pseudogene build log
│   ├── 02_genome_download.log          # Genome download log
│   ├── validation_checklist.log        # Validation results ✅
│   └── [...environment logs...]
│
└── outputs/                            # Test outputs
```

---

## Key Findings

### Bowtie 1 + Bowtie 2 Compatibility ✅
- **Discovery**: Bowtie 1 can read Bowtie 2 `.bt2` index files directly
- **Benefit**: Use pre-built Bowtie 2 indices (fast download) instead of building from FASTA (1-3 hours)
- **Tested**: All 3 genome indices tested and working with Bowtie 1 queries
- **Download time**: 12 minutes total (vs 1-3 hours for building from scratch)

### File Formats
- **Pseudogene indices**: Bowtie 1 (`.ebwt`) — built from shipped FASTA files
- **Genome indices**: Bowtie 2 (`.bt2`) — downloaded pre-built from AWS

Both formats work seamlessly with the ProbeDesign pipeline.

---

## Troubleshooting

### "command not found: probedesign"
```bash
micromamba activate probedesign
```

### "Index not found"
```bash
ls -1 bowtie_indexes/*.ebwt | wc -l  # Should show: 18 (pseudogene)
ls -1 bowtie_indexes/*.bt2 | wc -l   # Should show: 18 (genome)
```

### Test indices
```bash
./jyl/setup/03_test_indices.sh
# Should show 6/6 indices PASS
```

---

## Documentation

| Document | Description |
|----------|-------------|
| [`docs/00_index_setup_guide.md`](docs/00_index_setup_guide.md) | Complete index setup guide |
| [`docs/01_environment_setup.md`](docs/01_environment_setup.md) | Environment setup details |
| [`../jyl_implementation_plan.md`](../jyl_implementation_plan.md) | Overall implementation plan |
| [`../jyl_design_parameters.md`](../jyl_design_parameters.md) | Pipeline parameters reference |
| [`../../CLAUDE.md`](../../CLAUDE.md) | Main project documentation |

---

## Next Steps

1. **Run validation tests**: `./run_tests.sh` (from repository root)
2. **Design probes**: Use `probedesign design` with your FASTA files
3. **Enable masking features**: Add `--pseudogene-mask` and `--genome-mask` flags

All setup complete! 🎉
