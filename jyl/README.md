# jyl Setup — Complete ✅

**Location**: `/Users/jefflee/Documents/Github/smFISHProbeDesign/jyl/`
**Last Updated**: 2026-02-17

---

## Status: All Complete ✅

| Component | Status | Time | Size |
|-----------|--------|------|------|
| Environment | ✅ Complete | 5 min | - |
| Pseudogene indices | ✅ Complete | 47 sec | 84 MB |
| Genome indices | ✅ Complete | 12 min | 6.7 GB |
| Testing | ✅ Passed | <1 min | - |

**Total setup time**: ~18 minutes
**Total disk usage**: ~6.8 GB (indices only)

---

## Quick Start

📖 **Main Guide**: [`setup/README.md`](setup/README.md) — Complete setup documentation

### Verify Setup:
```bash
micromamba activate probedesign
probedesign --version  # Should show: 0.1.0
bowtie --version       # Should show: bowtie-align-s version 1.3.1
```

### Test Indices:
```bash
./jyl/setup/03_test_indices.sh
# Should show: All 6 indices PASS
```

### Run Validation:
```bash
./run_tests.sh
# All tests should pass
```

---

## Directory Structure

```
jyl/
├── README.md                         # This file
├── jyl_design_parameters.md          # Pipeline parameters reference
├── jyl_implementation_plan.md        # Overall implementation plan
│
└── setup/                            # All setup scripts & docs ✅
    ├── README.md                     # Main setup guide
    │
    ├── 01_build_pseudogene_indices.sh  # Build pseudogene indices ✅
    ├── 02_download_genome_indices.sh   # Download genome indices ✅
    ├── 03_test_indices.sh              # Test all indices ✅
    │
    ├── docs/                         # Documentation
    │   ├── 00_index_setup_guide.md   # Complete index guide
    │   ├── 01_environment_setup.md   # Environment setup
    │   └── 01_environment_setup_log.md
    │
    └── logs/                         # Build logs
        ├── 01_pseudogene_build.log   # Pseudogene build ✅
        ├── 02_genome_download.log    # Genome download ✅
        └── [...environment logs...]
```

---

## What's Ready

### ✅ All Features Available

- **Basic probe design**: `probedesign design input.fa --probes 32`
- **Repeat masking**: `--repeatmask-file masked.fa`
- **Pseudogene masking**: `--pseudogene-mask`
- **Genome masking**: `--genome-mask`

### ✅ All Indices Built and Tested

**Pseudogene Indices** (Bowtie 1 `.ebwt`):
- Human: 18,046 sequences, 35 MB
- Mouse: 19,086 sequences, 39 MB
- Drosophila: 2,204 sequences, 10 MB

**Genome Indices** (Bowtie 2 `.bt2`, compatible with Bowtie 1):
- Human GRCh38: 3.5 GB
- Mouse mm10: 3.1 GB
- Drosophila BDGP6: 176 MB

All 6 indices tested and working ✅

---

## Key Finding: Bowtie 1 ↔ Bowtie 2 Compatibility

**Discovery**: Bowtie 1 can read Bowtie 2 `.bt2` index files directly!

**Impact**:
- Download pre-built Bowtie 2 indices (12 minutes) instead of building from FASTA (1-3 hours)
- Saves significant setup time
- All genome indices use this approach

**Tested**: All 3 genome indices work perfectly with Bowtie 1 queries ✅

---

## Usage Examples

### Basic probe design:
```bash
probedesign design input.fa --probes 32
```

### With pseudogene masking:
```bash
probedesign design input.fa --probes 32 --pseudogene-mask
```

### With all masking features:
```bash
probedesign design input.fa --probes 32 \
  --repeatmask-file masked.fa \
  --pseudogene-mask \
  --genome-mask
```

---

## Validation

Run full pipeline tests:
```bash
./run_tests.sh
```

**Expected Results**:
- Test 1 (CDKN1A_32): 100% match ✅
- Test 2 (KRT19): 100% match ✅
- Test 3 (EIF1_HCR): 78% match ✅

---

## Documentation

| Document | Description |
|----------|-------------|
| [`setup/README.md`](setup/README.md) | **START HERE** — Complete setup guide |
| [`setup/docs/00_index_setup_guide.md`](setup/docs/00_index_setup_guide.md) | Detailed index setup documentation |
| [`jyl_implementation_plan.md`](jyl_implementation_plan.md) | Overall implementation plan |
| [`jyl_design_parameters.md`](jyl_design_parameters.md) | Pipeline parameters reference |
| [`../CLAUDE.md`](../CLAUDE.md) | Main project documentation |

---

## Troubleshooting

### "command not found: probedesign"
```bash
micromamba activate probedesign
```

### Test indices
```bash
./jyl/setup/03_test_indices.sh
# Should show: 6/6 indices PASS
```

### Check disk usage
```bash
du -sh bowtie_indexes/
# Should show: ~6.8 GB
```

---

**All setup complete and tested!** 🎉

See [`setup/README.md`](setup/README.md) for detailed documentation.
