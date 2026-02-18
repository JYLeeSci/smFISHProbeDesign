# ProbeDesign Setup Summary

**Date**: 2026-02-17
**Setup Duration**: Complete
**Platform**: macOS (Darwin 23.5.0)

## Overview

Successfully completed environment setup, package installation, bowtie verification, and minimal test (6.1) for the ProbeDesign Python pipeline. The Python implementation has been verified to produce **identical results** (100% match) to the MATLAB reference implementation.

## What Was Accomplished

### 1. Directory Structure ✓

Created organized setup folder:
```
jyl/setup/
├── scripts/          # 4 markdown documents with commands and explanations
│   ├── 01_create_environment.md
│   ├── 02_install_package.md
│   ├── 03_verify_installation.md
│   └── 04_minimal_test.md
├── logs/             # 6 command output logs
│   ├── 01_create_environment.log
│   ├── 02_install_package.log
│   ├── 03_verify_probedesign.log
│   ├── 04_verify_bowtie_path.log
│   ├── 05_verify_bowtie_version.log
│   └── 06_minimal_test.log
├── outputs/          # Test results
│   ├── CDKN1A_test_oligos.txt
│   ├── CDKN1A_test_seq.txt
│   ├── actual.txt
│   └── expected.txt
├── SETUP_LOG.md      # Master log with all details
└── README.md         # This summary
```

### 2. Micromamba Environment ✓

Created `probedesign` environment with:

| Package | Version | Source | Purpose |
|---------|---------|--------|---------|
| Python | 3.11.14 | conda-forge | Core language |
| click | 8.3.1 | conda-forge | CLI framework |
| bowtie | 1.3.1 | bioconda | Sequence alignment |
| pytest | 9.0.2 | conda-forge | Testing framework |

**Environment location**: `/Users/jefflee/micromamba/envs/probedesign`

**Total packages**: 37 (including dependencies)

### 3. ProbeDesign Package ✓

Installed probedesign v0.1.0 in editable mode:
- Installation method: `pip install -e .`
- Entry point: `probedesign` command
- Package location: `src/probedesign/`

### 4. Verification ✓

All components verified working:

| Component | Expected | Actual | Status |
|-----------|----------|--------|--------|
| probedesign version | 0.1.0 | 0.1.0 | ✓ PASS |
| bowtie path | In probedesign env | `/Users/jefflee/micromamba/envs/probedesign/bin/bowtie` | ✓ PASS |
| bowtie version | 1.3.1 (align-s) | 1.3.1 (bowtie-align-s, 64-bit) | ✓ PASS |

### 5. Minimal Test (6.1) ✓

**Test**: CDKN1A_32 with repeat masking from file

**Command**:
```bash
probedesign design test_cases/CDKN1A_32/CDKN1A.fa \
  --probes 32 \
  --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa \
  -o jyl/setup/outputs/CDKN1A_test
```

**Results**:

| Metric | Expected | Actual | Status |
|--------|----------|--------|--------|
| Probes generated | 32 | 32 | ✓ PASS |
| Sequence match | 100% | 100% | ✓ PASS |
| Repeat positions masked | 174 | 174 | ✓ PASS |
| Average badness score | ~0.08-0.09 | 0.0859 | ✓ PASS |

**Validation**: Python implementation produces **identical probes** to MATLAB reference implementation.

## How to Use

### Activate Environment

```bash
micromamba activate probedesign
```

### Run ProbeDesign

```bash
probedesign design input.fa --probes 32
```

### Run Without Activation

```bash
micromamba run -n probedesign probedesign design input.fa --probes 32
```

## What's Next

The following were intentionally **skipped for now** per your instructions:

### Reference Data & Indices (Section 5)

Not yet completed but documented in implementation plan:

1. **Pseudogene Indexes** (Section 5B) - ~3-5 min per species
   - Human: `bowtie-build probedesign/pseudogeneDBs/human.fasta bowtie_indexes/humanPseudo`
   - Mouse: `bowtie-build probedesign/pseudogeneDBs/mouse.fasta bowtie_indexes/mousePseudo`
   - Drosophila: `bowtie-build probedesign/pseudogeneDBs/drosophila.fasta bowtie_indexes/drosophilaPseudo`

2. **Genome Indexes** (Section 5C) - Large downloads
   - Human GRCh38: ~2.7 GB compressed, ~3 GB extracted
   - Mouse mm10: ~2.5 GB compressed, ~2.7 GB extracted
   - Drosophila dm6: Build from FASTA (~165 MB final)

3. **RepeatMasker & Dfam Databases** (Section 5D) - Optional, large
   - RepeatMasker installation: `micromamba install repeatmasker`
   - Dfam partition 7 (mammals): ~8.9 GB compressed, ~56 GB extracted
   - Dfam partition 1 (drosophila): ~7.8 GB compressed, ~45 GB extracted

### Additional Tests (Section 6.2-6.4)

Require bowtie indexes to be built first:
- **Test 6.2**: KRT19 with pseudogene + genome masking (expected: 100% match)
- **Test 6.3**: EIF1 HCR probes with 52bp oligos (expected: ~78% match)
- **Test 6.4**: Automated test script `./run_tests.sh` (runs all tests)

## Files Created

### Documentation (Markdown)
- `SETUP_LOG.md` - Comprehensive setup log with all steps
- `scripts/01_create_environment.md` - Environment creation details
- `scripts/02_install_package.md` - Package installation details
- `scripts/03_verify_installation.md` - Verification details
- `scripts/04_minimal_test.md` - Test execution and validation
- `README.md` - This summary

### Logs (Command Outputs)
- `logs/01_create_environment.log` - Micromamba environment creation output
- `logs/02_install_package.log` - pip install output
- `logs/03_verify_probedesign.log` - Version verification
- `logs/04_verify_bowtie_path.log` - Bowtie path check
- `logs/05_verify_bowtie_version.log` - Bowtie version output
- `logs/06_minimal_test.log` - Test execution output

### Test Outputs
- `outputs/CDKN1A_test_oligos.txt` - Generated 32 probes
- `outputs/CDKN1A_test_seq.txt` - Visual sequence alignment
- `outputs/actual.txt` - Sorted probe sequences (for diff)
- `outputs/expected.txt` - Reference probe sequences (for diff)

## Key Findings

1. ✓ **Micromamba is correctly installed** and functional
2. ✓ **Bioconda bowtie** (NOT homebrew) is correctly installed from bioconda
3. ✓ **Python implementation is correct** - 100% match with MATLAB on minimal test
4. ✓ **Core algorithms verified**:
   - FASTA parsing
   - Thermodynamic calculations (Sugimoto 1995 parameters)
   - Repeat masking from file
   - Dynamic programming optimization
   - Probe selection and output formatting

## Implementation Plan Status

Updated `jyl_implementation_plan.md` with:
- Setup status section showing completed and pending steps
- Link to setup documentation in `jyl/setup/`
- Clear checkboxes for tracking progress

## Questions or Issues?

If you encounter any issues or have questions about:
- Building bowtie indexes
- Installing RepeatMasker
- Running additional tests
- Any error messages

Please refer to:
1. `SETUP_LOG.md` - Detailed logs and troubleshooting
2. `jyl_implementation_plan.md` - Full installation plan with all sections
3. Root `CLAUDE.md` - Project context and architecture

## Quick Reference

### Environment Info
```bash
micromamba info -e                          # List environments
micromamba list -n probedesign              # List packages
micromamba activate probedesign             # Activate environment
```

### ProbeDesign Commands
```bash
probedesign --help                          # Show help
probedesign --version                       # Show version
probedesign design --help                   # Show design command help
```

### Test Commands
```bash
# Minimal test (no indexes needed)
probedesign design test_cases/CDKN1A_32/CDKN1A.fa \
  --probes 32 \
  --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa

# With bowtie masking (requires indexes)
probedesign design test_cases/KRT19_withUTRs/KRT19_withUTRs.fa \
  --probes 32 \
  --pseudogene-mask \
  --genome-mask \
  --index-dir bowtie_indexes
```

---

**Setup completed successfully**: 2026-02-17
