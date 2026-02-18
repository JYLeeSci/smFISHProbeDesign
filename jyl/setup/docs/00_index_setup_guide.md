# Bowtie Index Setup — Documentation

**Date**: 2026-02-17
**Location**: `/Users/jefflee/Documents/Github/smFISHProbeDesign/jyl/`

## Summary

This document describes the setup of bowtie1 indices required for the ProbeDesign pipeline's masking features (pseudogene and genome masking).

## Status

### ✅ Completed

- **Pseudogene indices** (3 species) ← **READY TO USE**
- Test scripts for index validation
- Documentation

### ⏳ Pending

- **Genome indices** (3 species) ← **REQUIRES BUILDING (1-3 hours)**

## What Was Done

### 1. Pseudogene Indices (✅ Built and Tested)

**Script**: `jyl/build_pseudogene_indices.sh`
**Log**: `jyl/pseudogene_index_build.log`

Built Bowtie 1 indices from shipped pseudogene FASTA files:

| Species | Index Prefix | Source FASTA | Sequences | Build Time | Index Size |
|---------|-------------|--------------|-----------|------------|-----------|
| Human | `humanPseudo` | `probedesign/pseudogeneDBs/human.fasta` | 18,046 | 23s | 35 MB |
| Mouse | `mousePseudo` | `probedesign/pseudogeneDBs/mouse.fasta` | 19,086 | 22s | 39 MB |
| Drosophila | `drosophilaPseudo` | `probedesign/pseudogeneDBs/drosophila.fasta` | 2,204 | 2s | 10 MB |

**Total**: 18 index files (6 per species), ~84 MB

**Index Files Created**:
```
bowtie_indexes/humanPseudo.{1,2,3,4,rev.1,rev.2}.ebwt
bowtie_indexes/mousePseudo.{1,2,3,4,rev.1,rev.2}.ebwt
bowtie_indexes/drosophilaPseudo.{1,2,3,4,rev.1,rev.2}.ebwt
```

**Verification**: All indices tested successfully with minimal bowtie queries (see `jyl/test_indices.sh`).

---

### 2. Genome Indices (⏳ Scripts Ready, Building Required)

**Script**: `jyl/build_genome_indices.sh`
**Expected Log**: `jyl/genome_index_build.log` (created when script runs)

#### Why Build from FASTA?

- Modern pre-built indices from AWS S3/JHU are **Bowtie 2** (`.bt2` files)
- ProbeDesign pipeline requires **Bowtie 1** (`.ebwt` files)
- FTP servers are unreliable/slow
- Building from FASTA ensures correct index format

#### To Build Genome Indices

The script downloads genome FASTAs from NCBI/Ensembl and builds Bowtie 1 indices:

```bash
cd /Users/jefflee/Documents/Github/smFISHProbeDesign
micromamba activate probedesign
./jyl/build_genome_indices.sh
```

**Expected Time**: 1-3 hours total

| Species | Index Prefix | Source | Download Time | Build Time | Final Size (est) |
|---------|-------------|--------|---------------|-----------|------------------|
| Human | `GCA_000001405.15_GRCh38_no_alt_analysis_set` | NCBI GRCh38 | 15-30 min | 45-90 min | ~3 GB |
| Mouse | `mm10` | Ensembl GRCm38 | 10-20 min | 45-90 min | ~2.7 GB |
| Drosophila | `drosophila` | Ensembl BDGP6.46 (dm6) | 2-5 min | 5-10 min | ~165 MB |

**Notes**:
- Human and mouse builds use `--threads 4` for faster indexing
- Drosophila builds quickly (~15 minutes total)
- FASTA files are deleted after index build to save space
- Script skips species if 6 `.ebwt` files already exist

**To run in background** (recommended for long builds):
```bash
nohup ./jyl/build_genome_indices.sh > jyl/genome_index_build_output.txt 2>&1 &
```

---

## Scripts Created

| Script | Purpose | Usage |
|--------|---------|-------|
| `build_pseudogene_indices.sh` | Build pseudogene indices from shipped FASTA files | Run once (already done) |
| `build_genome_indices.sh` | Download genome FASTAs and build Bowtie 1 indices | Run when ready (1-3 hours) |
| `test_indices.sh` | Test all indices with minimal bowtie queries | Run after each build |

All scripts:
- Check for conda environment activation
- Log detailed output to `jyl/*.log` files
- Skip species if indices already exist (safe to re-run)
- Support color-coded terminal output

---

## Testing

### Pseudogene Index Tests (✅ Passed)

**Script**: `jyl/test_indices.sh`

```bash
cd /Users/jefflee/Documents/Github/smFISHProbeDesign
micromamba activate probedesign
./jyl/test_indices.sh
```

**Results**:
- ✅ Human pseudogene: 5 hits for test sequence `ATCGATCGATCGATCG`
- ✅ Mouse pseudogene: 7 hits for test sequence
- ✅ Drosophila pseudogene: 5 hits for test sequence

All pseudogene indices respond correctly to bowtie queries.

### Genome Index Tests (⏳ Pending Build)

Run `./jyl/test_indices.sh` after building genome indices to verify they work.

---

## .gitignore Verification (✅)

Checked `.gitignore` at repository root:

```gitignore
# Bowtie indexes (large genome files)
bowtie_indexes/
*.ebwt
*.ebwtl
```

**Result**: ✅ All index files are excluded from git by pattern matching.

**Disk Usage Check**:
```bash
du -sh bowtie_indexes/
# Current: ~84 MB (pseudogene indices only)
# After genome build: ~6 GB (pseudogene + 3 genomes)
```

---

## Directory Structure

```
smFISHProbeDesign/
├── jyl/
│   ├── build_pseudogene_indices.sh       ← Built pseudogene indices
│   ├── build_genome_indices.sh           ← Build genome indices (ready to run)
│   ├── test_indices.sh                   ← Test all indices
│   ├── pseudogene_index_build.log        ← Build log (pseudogenes)
│   ├── genome_index_build.log            ← Build log (genomes, after running script)
│   ├── jyl_implementation_plan.md        ← Overall setup plan
│   ├── jyl_design_parameters.md          ← Pipeline parameter documentation
│   └── jyl_index_setup.md                ← This file
│
├── bowtie_indexes/                       ← Index directory (in .gitignore)
│   ├── humanPseudo.*.ebwt                ← 6 files, 35 MB
│   ├── mousePseudo.*.ebwt                ← 6 files, 39 MB
│   ├── drosophilaPseudo.*.ebwt           ← 6 files, 10 MB
│   ├── GCA_000001405.15_GRCh38_no_alt_analysis_set.*.ebwt  ← 6 files (after build), ~3 GB
│   ├── mm10.*.ebwt                       ← 6 files (after build), ~2.7 GB
│   └── drosophila.*.ebwt                 ← 6 files (after build), ~165 MB
│
└── probedesign/
    └── pseudogeneDBs/                    ← Source FASTA files (shipped with repo)
        ├── human.fasta                   ← 26 MB
        ├── mouse.fasta                   ← 28 MB
        └── drosophila.fasta              ← 2.1 MB
```

---

## Pipeline Integration

### How Indices Are Used

**Pseudogene masking** (`--pseudogene-mask`):
- Aligns 16-mers from input sequence against pseudogene index
- Bowtie command: `bowtie -f -v 0 -k 1 --quiet <pseudogene_index>`
- Masks positions with any hit (exact match)

**Genome masking** (`--genome-mask`):
- Aligns 12/14/16-mers against genome index at multiple stringencies
- Bowtie commands:
  - 12-mer: `bowtie -f -v 0 -k 5000` (threshold: 4000 hits)
  - 14-mer: `bowtie -f -v 0 -k 1000` (threshold: 500 hits)
  - 16-mer: `bowtie -f -v 0 -k 100` (threshold: 20 hits)
- Masks positions with hits above threshold

**Code References**:
- `src/probedesign/masking.py:251-272` — Pseudogene index mapping
- `src/probedesign/masking.py:297-326` — Genome index mapping
- `src/probedesign/masking.py:22-63` — Bowtie executable finding

---

## Next Steps

### To Complete Setup:

1. **Build genome indices** (when ready for 1-3 hour wait):
   ```bash
   cd /Users/jefflee/Documents/Github/smFISHProbeDesign
   micromamba activate probedesign
   ./jyl/build_genome_indices.sh
   ```

2. **Test genome indices**:
   ```bash
   ./jyl/test_indices.sh
   ```

3. **Run validation tests** with both masking types:
   ```bash
   # KRT19 test (requires pseudogene + genome indices)
   probedesign design test_cases/KRT19_withUTRs/KRT19_withUTRs.fa \
     --probes 32 --pseudogene-mask --genome-mask --index-dir bowtie_indexes \
     -o /tmp/KRT19_test

   # Compare with reference
   awk '{print $5}' /tmp/KRT19_test_oligos.txt | sort > /tmp/actual.txt
   awk '{print $5}' test_cases/KRT19_withUTRs/KRT19_withUTRs_oligos.txt | sort > /tmp/expected.txt
   diff /tmp/actual.txt /tmp/expected.txt
   # Expect: no output (100% match)
   ```

4. **Run automated test suite**:
   ```bash
   ./run_tests.sh
   # All tests should pass after genome indices are built
   ```

---

## Troubleshooting

### "bowtie-build not found"
```bash
micromamba activate probedesign
which bowtie-build  # Should show: /Users/jefflee/micromamba/envs/probedesign/bin/bowtie-build
```

### "SKIP: Index not found"
Check that 6 `.ebwt` files exist for the index:
```bash
ls -1 bowtie_indexes/<index_prefix>.*.ebwt | wc -l
# Should return: 6
```

### Disk space issues
```bash
df -h .
# Need ~7 GB free for all indices (84 MB pseudogene + ~6 GB genome)
```

### Slow genome build
Genome index building is CPU-intensive. Expected times:
- Drosophila: 15 minutes total
- Mouse: 1-2 hours
- Human: 1.5-2.5 hours

Use `nohup ... &` to background the process.

---

## References

- **CLAUDE.md** — Main project documentation (see "Bowtie Setup" section)
- **jyl_implementation_plan.md** — Overall setup plan (Section 5: Reference Data & Indices)
- **BOWTIE.md** — Detailed bowtie installation and index setup instructions
- **run_tests.sh** — Automated test suite for validation

---

## Log Files

All build and test operations log to files in `jyl/`:

- `pseudogene_index_build.log` — Pseudogene index builds (✅ complete)
- `genome_index_build.log` — Genome index builds (created when script runs)
- Output from test scripts is printed to terminal (not logged)

To review pseudogene build log:
```bash
cat jyl/pseudogene_index_build.log
```

---

**Last Updated**: 2026-02-17
**Status**: Pseudogene indices ready ✅ | Genome indices pending ⏳ (script ready to run)
