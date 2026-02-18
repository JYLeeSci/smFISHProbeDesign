# ProbeDesign Setup Log

**Date**: 2026-02-17
**Platform**: macOS (Darwin 23.5.0)
**Micromamba**: `/Users/jefflee/.local/bin/micromamba`

## Setup Progress

### 1. Setup Folder Structure ✓

Created the following directory structure:
```
jyl/setup/
├── scripts/      # Setup scripts and commands
├── logs/         # Command outputs and logs
└── outputs/      # Test outputs and results
```

### 2. Environment Setup ✓

Created micromamba environment named `probedesign` with:
- Python 3.11.14
- click 8.3.1
- bowtie 1.3.1
- pytest 9.0.2

**Status**: Complete

**Environment location**: `/Users/jefflee/micromamba/envs/probedesign`

**Log**: `logs/01_create_environment.log`

**Script**: `scripts/01_create_environment.md`

### 3. Package Installation ✓

Installed probedesign package in editable mode:
```bash
pip install -e .
```

**Version**: 0.1.0

**Status**: Complete

**Log**: `logs/02_install_package.log`

**Script**: `scripts/02_install_package.md`

### 4. Verification ✓

Verified all installations:

| Component | Expected | Actual | Status |
|-----------|----------|--------|--------|
| probedesign version | 0.1.0 | 0.1.0 | ✓ PASS |
| bowtie path | In env | `/Users/jefflee/micromamba/envs/probedesign/bin/bowtie` | ✓ PASS |
| bowtie version | 1.3.1 | 1.3.1 (bowtie-align-s) | ✓ PASS |

**Logs**:
- `logs/03_verify_probedesign.log`
- `logs/04_verify_bowtie_path.log`
- `logs/05_verify_bowtie_version.log`

**Script**: `scripts/03_verify_installation.md`

### 5. Minimal Test (6.1) ✓

Ran CDKN1A test without bowtie indexes:
```bash
probedesign design test_cases/CDKN1A_32/CDKN1A.fa \
  --probes 32 \
  --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa \
  -o jyl/setup/outputs/CDKN1A_test
```

**Results**:
- Probes generated: 32/32
- Match with reference: 100%
- Repeat positions masked: 174
- Average badness score: 0.0859

**Status**: ✓ PASS - Python implementation matches MATLAB reference exactly

**Output files**:
- `outputs/CDKN1A_test_oligos.txt`
- `outputs/CDKN1A_test_seq.txt`

**Log**: `logs/06_minimal_test.log`

**Script**: `scripts/04_minimal_test.md`

## Summary

All setup and validation steps completed successfully:

- [x] Folder structure created
- [x] Micromamba environment created
- [x] ProbeDesign package installed
- [x] Bowtie verified (correct version from bioconda)
- [x] Minimal test passed (100% match)

## Next Steps (Not Yet Completed)

### Reference Data & Indices (Section 5)

The following are needed for full functionality but were **skipped for now**:

1. **Pseudogene Indexes** (Section 5B)
   - Human: `bowtie-build probedesign/pseudogeneDBs/human.fasta bowtie_indexes/humanPseudo`
   - Mouse: `bowtie-build probedesign/pseudogeneDBs/mouse.fasta bowtie_indexes/mousePseudo`
   - Drosophila: `bowtie-build probedesign/pseudogeneDBs/drosophila.fasta bowtie_indexes/drosophilaPseudo`

2. **Genome Indexes** (Section 5C)
   - Human GRCh38: ~2.7 GB download
   - Mouse mm10: ~2.5 GB download
   - Drosophila dm6: Build from FASTA

3. **RepeatMasker & Dfam Databases** (Section 5D)
   - Install RepeatMasker (optional, for `--repeatmask` flag)
   - Download species-specific Dfam partitions (~56 GB for mammals)

### Additional Tests

Once indexes are built:
- Test 6.2: KRT19 with pseudogene + genome masking
- Test 6.3: EIF1 HCR probes
- Test 6.4: Automated test script (`./run_tests.sh`)

## Activation Instructions

To use the installed environment:

```bash
micromamba activate probedesign
probedesign design input.fa --probes 32
```

Or run commands without activating:

```bash
micromamba run -n probedesign probedesign design input.fa --probes 32
```
