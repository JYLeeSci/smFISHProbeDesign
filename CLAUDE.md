# Claude Code Context for ProbeDesign

This file provides context for Claude Code when working on this project.

## Project Overview

ProbeDesign is a tool for designing oligonucleotide probes for single molecule RNA FISH experiments. It was originally written in MATLAB and has been ported to Python for easier installation and use.

**Key goal**: The Python implementation should produce **identical results** to the MATLAB reference implementation (`probedesign/findprobesLocal.m`).

For python implementation remember to use the user's virtual environment and search for micromamba, mamba or conda probedesign environment, and avoid installing python packages to user's base CLI environment. 

## Architecture

### Python Package (`src/probedesign/`)

| File | Purpose |
|------|---------|
| `cli.py` | Click-based command-line interface |
| `core.py` | Main probe design algorithm (badness calculation, DP optimization, mixed-length support) |
| `thermodynamics.py` | Gibbs free energy and melting temperature calculations (Sugimoto 1995 params) |
| `masking.py` | Bowtie-based sequence masking (pseudogene, genome), RepeatMasker integration, mixed-length mask support |
| `sequence.py` | Sequence utilities (reverse complement, GC%, validation) |
| `fasta.py` | FASTA file I/O with junction marker handling |
| `output.py` | Output file generation (_oligos.txt, _seq.txt) |

### Streamlit App (`streamlit_app/`)

| File | Purpose |
|------|---------|
| `app.py` | Main Streamlit UI — sidebar parameters, single/batch mode panels, result display |
| `utils.py` | Backend helpers — FASTA validation, `run_design()` wrapper, `run_batch()` runner, prerequisite checks |
| `README.md` | End-user installation and usage guide |

**Running the app**:
```bash
cd streamlit_app
streamlit run app.py
```

**Features**:
- **Single sequence mode**: Upload or paste a FASTA file, design probes interactively
- **Batch mode**: Process multiple FASTA files from a directory or upload, write results to an output directory
- **Mixed-length support**: Toggle "Mixed range" in sidebar to use variable-length probes (e.g., 18-22bp)
- **Masking options**: Pseudogene mask, genome mask, RepeatMasker (auto or file)
- **Downloads**: `_oligos.txt`, `_seq.txt`, bowtie hit files, batch summary TSV

**Parameter flow**:
```
Sidebar UI (app.py) → run_design() (utils.py) → design_probes() (core.py)
Sidebar UI (app.py) → run_batch() (utils.py) → run_design() per file → design_probes()
```

**Key functions in `utils.py`**:

| Function | Purpose |
|----------|---------|
| `run_design()` | Wraps `design_probes()` with stdout capture, RepeatMasker auto-execution, error handling |
| `run_batch()` | Iterates FASTA files, calls `run_design()` per file, writes output files, tracks progress |
| `check_prerequisites()` | Verifies bowtie, RepeatMasker, and index files are available |
| `validate_fasta_text()` | Validates pasted FASTA format (headers, valid bases) |
| `clean_name_for_prefix()` | Converts filenames/headers to clean output name prefixes |

### MATLAB Code (`probedesign/`)

| File | Purpose |
|------|---------|
| `findprobesLocal.m` | Main MATLAB function - **reference implementation** |
| `thermo_RNA_DNA.m` | Thermodynamics (Sugimoto 1995 parameters) |
| `find_best_matches.m` | Dynamic programming algorithm |
| `pseudogeneDBs/` | Pseudogene FASTA files for masking |

## Algorithm

1. **Read input**: Parse FASTA, concatenate multi-entry files with `>` junction markers
2. **Calculate badness**: For each position, compute `(gibbs - target)^2`. Positions with invalid chars or out-of-range Gibbs get `inf` badness. In mixed-length mode (`-l 18-22`), badness is calculated for each (position, length) pair as a 2D table.
3. **Apply masking**:
   - R mask: Repeat regions (from N's in input or separate repeatmask file)
   - P mask: Pseudogene alignments (bowtie to pseudogene DBs)
   - B mask: Genome alignments (bowtie to genome with multiple stringencies)
   - F mask: Thermodynamic filtering (badness == inf)
4. **Dynamic programming**: Find optimal probe placements minimizing average badness. Fixed-length mode uses start-position DP; mixed-length mode uses end-position DP to handle variable probe spacing.
5. **Generate output**: Create oligo and sequence alignment files

## Output Format

### `_oligos.txt`
Tab-separated: `index  GC%  Tm  Gibbs  sequence  name`

### `_seq.txt`
Visual alignment showing:
- Line 1: Original sequence
- Line 2: R mask (if repeat masking enabled)
- Line 3+: P mask, B mask (if enabled)
- Last mask line: F mask (thermodynamic filtering)
- Probe annotations below each block

## Test Cases (`test_cases/`)

| Test Case | Description | Expected Match |
|-----------|-------------|----------------|
| `CDKN1A_32/` | 32 probes, manual repeat masking (`--repeatmask-file`) | 100% |
| `CDKN1A_32/` | 32 probes, automatic RepeatMasker (`--repeatmask`) | 100% |
| `KRT19_withUTRs/` | 6 probes, pseudogene + genome masking | 100% |
| `EIF1_CDS_HCR/` | HCR probes (52bp oligos) | ~78% (partial expected) |
| `mixed_length_test/` | Synthetic GC-variable sequence for mixed-length testing | N/A (new feature) |

### Running Tests

Use the automated test script:

```bash
./run_tests.sh
```

This will run all tests and report pass/fail. Tests require bowtie installed via conda (NOT homebrew).

### Manual Test Commands

```bash
# Test 1: CDKN1A with repeatmask-file (100% match expected)
probedesign design test_cases/CDKN1A_32/CDKN1A.fa --probes 32 \
  --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa

# Test 2: CDKN1A with automatic RepeatMasker (100% match expected, requires partition 7)
probedesign design test_cases/CDKN1A_32/CDKN1A.fa --probes 32 --repeatmask

# Test 3: KRT19 with bowtie masking (100% match expected)
probedesign design test_cases/KRT19_withUTRs/KRT19_withUTRs.fa --probes 32 \
  --pseudogene-mask --genome-mask --index-dir bowtie_indexes

# Test 4: EIF1 HCR probes (78% match expected)
probedesign design test_cases/EIF1_CDS_HCR/EIF1_Exons.fasta --probes 20 \
  -l 52 --target-gibbs -60 --allowable-gibbs -80,-40 \
  --pseudogene-mask --genome-mask --index-dir bowtie_indexes

# Test 5: Mixed-length probes on CDKN1A
probedesign design test_cases/CDKN1A_32/CDKN1A.fa --probes 48 \
  --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa -l 18-22

# Test 6: Mixed-length on synthetic GC-variable sequence
probedesign design test_cases/mixed_length_test/mixed_gc.fa --probes 20 -l 18-22
```

## Important Implementation Details

### Thermodynamics

- Uses Sugimoto 1995 RNA-DNA hybrid parameters
- Salt concentration: 0.33 M
- Primer concentration: 50 µM
- Default target Gibbs: -23 kcal/mol (range: -26 to -20)

### F Mask Logic

The F mask shows thermodynamic filtering - positions where a probe **cannot start**:
- `F` at position i if `badness[i] == inf` OR `i >= goodlen` (past valid positions)
- Sequence char at position i if `badness[i]` is finite (valid probe start)

This matches MATLAB `mask_string(inseq, badness==inf, 'F')`.

### Mixed-Length Probe Design

Allows probes of varying lengths (e.g., 18-22bp) instead of a single fixed length. This increases probe count on sequences with variable GC content, where GC-rich regions need shorter probes and AT-rich regions need longer probes to meet the Gibbs FE target.

**CLI usage**:
```bash
# Fixed length (default, backward compatible)
probedesign design input.fa -l 20

# Mixed length range
probedesign design input.fa -l 18-22
```

**Streamlit app**: Select "Mixed range" under "Oligo length mode" in the sidebar, then set min/max lengths (e.g., 18 and 22).

**How it works**:
1. `calculate_badness_mixed()` builds a 2D table `badness[pos][len_idx]` for all positions and all candidate lengths
2. `find_best_probes_mixed()` uses an end-position DP (`dp[e][k]` = best score for k+1 probes ending at or before position e) to handle variable spacing between probes of different lengths
3. The DP considers ALL valid lengths at every position and finds the globally optimal combination
4. The same `target_gibbs` and `allowable_gibbs` are used for all lengths (ensures similar binding thermodynamics)

**Key functions**:

| Function | File | Purpose |
|----------|------|---------|
| `calculate_badness_mixed()` | `core.py` | 2D badness table for all (position, length) pairs |
| `find_best_probes_mixed()` | `core.py` | End-position DP for variable-length probe placement |
| `mask_to_badness_mixed()` | `masking.py` | Convert position-level mask to 2D probe-level badness |

**F mask in mixed-length mode**: A position shows 'F' only if NO length in the range produces a valid probe (all lengths give `inf` badness).

**Backward compatibility**: When `-l 20` (single integer) is used, the original fixed-length code path runs unchanged. The mixed-length code is only activated by range syntax (`-l 18-22`).

### Junction Handling

Multi-entry FASTA files use `>` to mark exon junctions. Probes cannot span junctions.

### Repeat Masking

Three modes:
1. **Auto-detect**: N's in input file are treated as masked regions
2. **Automatic**: `--repeatmask` runs RepeatMasker locally (requires installation)
3. **Manual file**: `--repeatmask-file` provides separate file with N's marking repeats

### RepeatMasker Setup

See [REPEATMASKER.md](REPEATMASKER.md) for installation instructions.

Key points:
- Install via **conda/mamba**: `mamba install -c bioconda -c conda-forge repeatmasker`
- Requires Dfam database partitions by species:
  - Human/Mouse/Rat: Partition 7 (Mammalia) - ~8.9GB compressed, **~56GB extracted**
- Location: `/opt/homebrew/Caskroom/miniforge/base/share/RepeatMasker/Libraries/famdb/`

#### RepeatMasker Functions in `masking.py`

| Function | Purpose |
|----------|---------|
| `find_repeatmasker()` | Locates RepeatMasker executable |
| `run_repeatmasker(fasta_path, species)` | Runs RepeatMasker, returns path to `.masked` file |
| `repeatmasker_mask_to_sequence(seq, masked_path)` | Converts masked output to binary mask |

#### Quick Test

```bash
# Test RepeatMasker integration (requires partition 7 installed)
probedesign design test_cases/CDKN1A_32/CDKN1A.fa --repeatmask --probes 32

# Compare with manual repeatmask-file (should produce identical results)
probedesign design test_cases/CDKN1A_32/CDKN1A.fa \
  --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa --probes 32
```

## Bowtie Setup

See [BOWTIE.md](BOWTIE.md) for full installation instructions.

Key points:
- Uses **Bowtie 1** (not Bowtie 2) for short read alignment
- Install via **conda/mamba** (NOT homebrew - that's a different tool)
- The test script auto-detects bowtie at `/opt/homebrew/Caskroom/miniforge/base/bin/bowtie`

### Required Indexes

Two types of indexes are needed in `bowtie_indexes/`:

1. **Pseudogene indexes** - Build from FASTA files in `probedesign/pseudogeneDBs/`:
   ```bash
   cd bowtie_indexes
   bowtie-build ../probedesign/pseudogeneDBs/human.fasta humanPseudo
   ```

2. **Genome indexes** - Download pre-built:
   ```bash
   cd bowtie_indexes
   curl -O ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/GRCh38_no_alt.zip
   unzip GRCh38_no_alt.zip
   ```

### Index Naming Convention

| Species | Pseudogene Index | Genome Index |
|---------|------------------|--------------|
| human | `humanPseudo` | `GCA_000001405.15_GRCh38_no_alt_analysis_set` |
| mouse | `mousePseudo` | `mm10` |

Verify bowtie installation:
```bash
bowtie --version  # Should show "bowtie-align-s version 1.x.x"
```

## Common Development Tasks

### Adding a new CLI option

1. Add Click option in `cli.py`
2. Pass parameter to `design_probes()` in `core.py`
3. Update `design_probes()` function signature and logic
4. Update README.md and this file

### Validating against MATLAB

1. Run MATLAB command with specific options
2. Save `_oligos.txt` and `_seq.txt` output
3. Run Python with equivalent options
4. Compare outputs with `diff`

### Debugging probe differences

1. Check F mask line matches (thermodynamic filtering)
2. Check mask strings (R, P, B) match
3. Verify badness calculation for specific positions
4. Compare probe positions and sequences

## Code Style

- Python 3.8+ with type hints
- Click for CLI
- No external dependencies except numpy (optional)
- Match MATLAB output format exactly for validation
