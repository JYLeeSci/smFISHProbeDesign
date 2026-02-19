# Claude Code Context for ProbeDesign

This file provides context for Claude Code when working on this project.

## Project Overview

ProbeDesign designs oligonucleotide probes for single molecule RNA FISH experiments. Originally written in MATLAB (`probedesign/findprobesLocal.m`), the active development is now the **Python implementation** with a **Streamlit web app** as the primary user interface.

**Development priorities**: Adding features to the Python implementation and Streamlit app. MATLAB parity is a secondary concern.

## Environment Setup

**Always use the user's virtual environment.** Search for `micromamba`, `mamba`, or `conda` probedesign environment. Never install Python packages to the user's base CLI environment.

```bash
# Activate environment
micromamba activate probedesign   # or: mamba / conda activate probedesign

# Setup from scratch
chmod +x setup_all.sh && ./setup_all.sh

# Launch web app
streamlit run streamlit_app/app.py

# CLI usage
probedesign design input.fa --probes 48
```

The `setup_all.sh` script (in repo root) auto-detects micromamba/mamba/conda, creates the environment from `environment.yml`, installs the package, builds pseudogene indices, and downloads genome indices.

## Repository Structure

```
smFISHProbeDesign/
├── src/probedesign/          # Python package (active development)
├── streamlit_app/            # Streamlit web app
├── probedesign/              # Legacy MATLAB code + pseudogeneDBs/
├── docs/                     # Documentation
├── test_cases/               # Validation test cases
├── bowtie_indexes/           # Genome/pseudogene indices (gitignored)
├── setup_all.sh              # Automated setup script
├── environment.yml           # Conda environment specification
├── pyproject.toml            # Python package config
└── README.md                 # User-facing documentation
```

## Architecture

### Python Package (`src/probedesign/`)

| File | Purpose |
|------|---------|
| `cli.py` | Click-based CLI; parses `--oligo-length` as int or range (e.g. `18-22`) |
| `core.py` | Main algorithm: badness calculation, DP optimization, mixed-length support |
| `thermodynamics.py` | Gibbs free energy and Tm calculations (Sugimoto 1995 params) |
| `masking.py` | Bowtie-based masking (pseudogene, genome), RepeatMasker, mixed-length mask support |
| `sequence.py` | Sequence utilities (reverse complement, GC%, validation) |
| `fasta.py` | FASTA file I/O with junction marker handling |
| `output.py` | Output file generation (_oligos.txt, _seq.txt) |

### Streamlit App (`streamlit_app/`)

| File | Purpose |
|------|---------|
| `app.py` | Main Streamlit UI — sidebar parameters, single/batch mode, result display |
| `utils.py` | Backend helpers — FASTA validation, `run_design()`, `run_batch()`, prerequisite checks |
| `README.md` | Quick launch instructions (full docs in root README) |

**Running the app**:
```bash
streamlit run streamlit_app/app.py
```

**Features**:
- **Single sequence mode**: Upload or paste a FASTA file, design probes interactively
- **Batch mode**: Process multiple FASTA files from a directory or upload, write results to an output directory
- **Mixed-length support**: Toggle "Mixed range" in sidebar to use variable-length probes (e.g. 18-22bp)
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
| `run_design()` | Wraps `design_probes()` with stdout capture, RepeatMasker auto-execution, error handling; accepts `mixed_lengths` parameter |
| `run_batch()` | Iterates FASTA files, calls `run_design()` per file, writes output files, tracks progress; extracts `mixed_lengths` from params dict |
| `check_prerequisites()` | Verifies bowtie, RepeatMasker, and index files are available |
| `validate_fasta_text()` | Validates pasted FASTA format (headers, valid bases) |
| `clean_name_for_prefix()` | Converts filenames/headers to clean output name prefixes |

### Documentation (`docs/`)

| File | Purpose |
|------|---------|
| `BOWTIE.md` | Bowtie installation and index setup |
| `REPEATMASKER.md` | RepeatMasker installation (~56 GB Dfam database) |
| `MIXED_LENGTH_PLAN.md` | Mixed-length probe design plan and analysis |
| `probe_design_principles.md` | Detailed algorithm reference (thermodynamics, masking, DP) |
| `summary.md` | Development notes (Bowtie 1/2 compatibility, installation findings) |

### Legacy Code

| Directory | Purpose |
|-----------|---------|
| `probedesign/` | MATLAB code (`findprobesLocal.m`) and pseudogene FASTA files (`pseudogeneDBs/`) |
| `DesignServer/` | Original Python server code (not actively maintained) |
| `maskprobes/` | Mask probe utilities (legacy) |
| `panprobedesign/` | Pan-probe MATLAB code (legacy) |

## Algorithm

1. **Read input**: Parse FASTA, concatenate multi-entry files with `>` junction markers
2. **Calculate badness**: For each position, compute `(gibbs - target)^2`. Positions with invalid chars or out-of-range Gibbs get `inf`. In mixed-length mode, computed for each (position, length) pair as a 2D table.
3. **Apply masking**:
   - R mask: Repeat regions (N's in input, separate file, or RepeatMasker)
   - P mask: Pseudogene alignments (bowtie to pseudogene DBs)
   - B mask: Genome alignments (bowtie to genome, tiered stringencies)
   - F mask: Thermodynamic filtering (badness == inf)
4. **Dynamic programming**: Find optimal probe placements minimizing average badness
   - Fixed-length: start-position DP with fixed spacing
   - Mixed-length: end-position DP to handle variable probe spacing
5. **Generate output**: Create `_oligos.txt` and `_seq.txt` files

## Mixed-Length Probe Design

Probes of varying lengths (e.g. 18-22bp) instead of a single fixed length. Maximises probe count on sequences with variable GC content — GC-rich regions get shorter probes, AT-rich regions get longer probes, all within the same Gibbs FE target range.

**CLI**: `probedesign design input.fa -l 18-22`

**Streamlit**: Select "Mixed range" under "Oligo length mode", set min/max lengths.

**How it works**:
1. `calculate_badness_mixed()` builds 2D table `badness[pos][len_idx]` for all positions × all lengths
2. `find_best_probes_mixed()` uses end-position DP: `dp[e][k]` = best score for k+1 probes ending at/before position e
3. Global DP considers ALL valid lengths at every position, finds optimal combination
4. Same `target_gibbs` and `allowable_gibbs` for all lengths

**Key functions**:

| Function | File | Purpose |
|----------|------|---------|
| `calculate_badness_mixed()` | `core.py` | 2D badness table |
| `find_best_probes_mixed()` | `core.py` | End-position DP |
| `mask_to_badness_mixed()` | `masking.py` | Position mask → 2D probe-level badness |

**F mask in mixed mode**: Position shows 'F' only if NO length gives finite badness.

**Backward compatibility**: `-l 20` (single int) runs the original fixed-length code path unchanged.

## Masking Details

### Bowtie Setup

- Uses **Bowtie 1** (not 2) for short read alignment — install via conda/bioconda, NOT homebrew
- Bowtie 1 v1.3+ reads Bowtie 2 `.bt2` index files directly (key finding — avoids multi-hour builds)
- `setup_all.sh` handles all index building/downloading automatically

### Required Indexes (in `bowtie_indexes/`)

| Type | Species | Index Name |
|------|---------|------------|
| Pseudogene | human | `humanPseudo` |
| Pseudogene | mouse | `mousePseudo` |
| Pseudogene | drosophila | `drosophilaPseudo` |
| Genome | human | `GCA_000001405.15_GRCh38_no_alt_analysis_set` |
| Genome | mouse | `mm10` |
| Genome | drosophila | `drosophila` |

### RepeatMasker

- Install via conda: `mamba install -c bioconda -c conda-forge repeatmasker`
- Requires Dfam database partitions: Human/Mouse/Rat need Partition 7 (~8.9GB compressed, ~56GB extracted)
- See [docs/REPEATMASKER.md](docs/REPEATMASKER.md)

## Output Format

### `_oligos.txt`
Tab-separated: `index  start  GC%  Tm  Gibbs  sequence  name`

### `_seq.txt`
Visual alignment: original sequence, R/P/B/F mask lines, probe annotations below each block.

## Test Cases (`test_cases/`)

| Test Case | Description | Expected |
|-----------|-------------|----------|
| `CDKN1A_32/` | 32 probes, manual repeat masking | 100% match |
| `KRT19_withUTRs/` | 6 probes, pseudogene + genome masking | 100% match |
| `EIF1_CDS_HCR/` | HCR probes (52bp oligos) | ~78% match |
| `mixed_length_test/` | Synthetic GC-variable, mixed-length | N/A (new feature) |

### Test Commands

```bash
# Fixed-length, repeat mask file
probedesign design test_cases/CDKN1A_32/CDKN1A.fa --probes 32 \
  --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa

# Bowtie masking
probedesign design test_cases/KRT19_withUTRs/KRT19_withUTRs.fa --probes 32 \
  --pseudogene-mask --genome-mask --index-dir bowtie_indexes

# Mixed-length
probedesign design test_cases/CDKN1A_32/CDKN1A.fa --probes 48 \
  --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa -l 18-22
```

## Thermodynamics

- Sugimoto 1995 RNA-DNA hybrid nearest-neighbor parameters
- Salt concentration: 0.33 M
- Primer concentration: 50 µM
- Default target Gibbs: -23 kcal/mol (range: -26 to -20)
- For HCR probes: `-l 52 --target-gibbs -60 --allowable-gibbs -80,-40`

## Common Development Tasks

### Adding a new CLI option

1. Add Click option in `cli.py`
2. Pass parameter to `design_probes()` in `core.py`
3. Update `design_probes()` function signature and logic
4. If applicable, add to Streamlit sidebar in `app.py` and pass through `utils.py`
5. Update README.md and this file

### Adding a new Streamlit feature

1. Add UI elements in `app.py` (sidebar or main panel)
2. Pass new parameter through `run_design()` in `utils.py`
3. Ensure batch mode also supports it via `run_batch()` params dict
4. Update README.md and this file

## Code Style

- Python 3.8+ with type hints
- Click for CLI, Streamlit for web app
- No external dependencies except numpy (optional) and streamlit/pandas (web app)
