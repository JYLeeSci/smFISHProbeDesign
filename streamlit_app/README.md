# ProbeDesign — User Guide

> Design oligonucleotide probes for single molecule RNA FISH experiments.
>
> This folder is the **single source of truth** for installation and usage.
> Start here regardless of your operating system or experience level.

---

## Table of Contents

1. [Before You Begin — Pick Your Setup Path](#1-before-you-begin--pick-your-setup-path)
2. [Step 1 — Install a Package Manager](#2-step-1--install-a-package-manager)
3. [Step 2 — Get the Repository](#3-step-2--get-the-repository)
4. [Step 3 — Run the Setup Script](#4-step-3--run-the-setup-script)
5. [Step 4 — Launch the Web App](#5-step-4--launch-the-web-app)
6. [Using the Command-Line Interface](#6-using-the-command-line-interface)
7. [Masking Options](#7-masking-options)
8. [Output Files](#8-output-files)
9. [Manual Environment Setup (Advanced)](#9-manual-environment-setup-advanced)
10. [Troubleshooting](#10-troubleshooting)
11. [Further Reading](#11-further-reading)

---

## 1. Before You Begin — Pick Your Setup Path

| Your OS | Supported? | Notes |
|---------|-----------|-------|
| **macOS** (Apple Silicon M1/M2/M3) | ✅ Full support | |
| **macOS** (Intel) | ✅ Full support | |
| **Linux** (x86_64 or aarch64) | ✅ Full support | |
| **Windows** | ⚠️ Via WSL2 only | [See Windows instructions below](#windows-users-wsl2-required) |

**What gets installed** — the setup script does all of this automatically:
- A conda environment (`probedesign`) with Python 3.11 and Bowtie 1.3.1
- The `probedesign` Python package + Streamlit web app
- ~84 MB of pseudogene alignment indices (built locally in ~1 min)
- ~6.7 GB of genome alignment indices (downloaded from AWS in ~12 min)

---

## 2. Step 1 — Install a Package Manager

You need **one** of: `micromamba`, `mamba`, or `conda`. The setup script auto-detects whichever you have.

### Option A — Miniforge (recommended for new users)

Miniforge installs both `conda` and `mamba` and is free for all uses.

1. Go to: **https://github.com/conda-forge/miniforge/releases/latest**
2. Download the installer for your system:

   | OS | File to download |
   |----|-----------------|
   | macOS Apple Silicon (M1/M2/M3) | `Miniforge3-MacOSX-arm64.sh` |
   | macOS Intel | `Miniforge3-MacOSX-x86_64.sh` |
   | Linux x86_64 | `Miniforge3-Linux-x86_64.sh` |
   | Linux ARM (e.g. Raspberry Pi) | `Miniforge3-Linux-aarch64.sh` |

3. Open **Terminal** (macOS/Linux) and run:
   ```bash
   bash ~/Downloads/Miniforge3-*.sh
   ```
   Follow the on-screen prompts. Accept the license and say **yes** when asked to initialize.

4. Close and reopen Terminal. Verify:
   ```bash
   conda --version
   # Expected: conda 24.x.x  (or similar)
   ```

### Option B — micromamba (faster, minimal install)

```bash
# macOS / Linux (one-liner):
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

Close and reopen Terminal. Verify:
```bash
micromamba --version
# Expected: 2.x.x
```

### Windows Users — WSL2 Required

Bowtie (the sequence alignment tool) is not available natively on Windows.
You must use **WSL2** (Windows Subsystem for Linux):

1. Open **PowerShell as Administrator** and run:
   ```powershell
   wsl --install
   ```
2. Restart your computer.
3. Open the newly installed **Ubuntu** app from the Start menu.
4. Inside Ubuntu, follow the **Linux** instructions above (Option A or B).
5. All subsequent commands in this guide should be run inside the Ubuntu terminal.

---

## 3. Step 2 — Get the Repository

If you haven't already cloned the repository:

```bash
# Install git if needed (skip if already installed):
#   macOS: git comes pre-installed; or: xcode-select --install
#   Linux: sudo apt install git   (Ubuntu/Debian)

git clone https://github.com/arjunrajlaboratory/ProbeDesign.git
cd ProbeDesign
```

If you downloaded a ZIP instead of using git:
1. Unzip it
2. Open Terminal and `cd` into the unzipped folder:
   ```bash
   cd ~/Downloads/ProbeDesign-main
   ```

---

## 4. Step 3 — Run the Setup Script

> **Run all commands from the repository root folder** (the folder containing `README.md`, `pyproject.toml`, etc.)

### Make the script executable (first time only)

```bash
chmod +x streamlit_app/setup_all.sh
```

### Run the full setup

```bash
./streamlit_app/setup_all.sh
```

The script will:
1. Auto-detect your package manager (micromamba / mamba / conda)
2. Create the `probedesign` environment from `streamlit_app/environment.yml`
3. Install the `probedesign` package and Streamlit
4. Build pseudogene alignment indices (~1 minute)
5. Download genome alignment indices (~12 minutes, ~6.7 GB)
6. Run 7 validation tests and report pass/fail

> **First run takes ~18 minutes** (mostly downloading genome indices). Subsequent runs are instant — the script skips anything already done.

### Skip genome downloads (optional)

If you only need pseudogene masking (not full genome masking), save ~12 minutes by skipping the 6.7 GB download:

```bash
./streamlit_app/setup_all.sh --skip-genome
```

### Expected output (success)

```
Detecting conda package manager…
Found: micromamba (2.x.x)

======================================================
  1/5  Create 'probedesign' environment
======================================================
...
======================================================
  5/5  Validation tests
======================================================

  [A] probedesign command works … PASS
  [B] bowtie command works … PASS
  [C] probe design — CDKN1A_32 (32 probes, 100% match expected) … PASS (32/32 match)
  [D] pseudogene index files present … PASS (18 files)
  [E] genome index files present … PASS (18 files)
  [F] bowtie queries against all 6 indices … PASS (6/6)
  [G] probe design — KRT19 (6 probes, pseudogene+genome mask, 100% match) … PASS (6/6 match)

  Tests: 7/7 passed
  All checks passed!
```

---

## 5. Step 4 — Launch the Web App

### Activate the environment

```bash
# Use whichever you installed:
micromamba activate probedesign
# OR:
mamba activate probedesign
# OR:
conda activate probedesign
```

Your terminal prompt will change to show `(probedesign)`.

### Launch Streamlit

```bash
streamlit run streamlit_app/app.py
```

Your browser opens automatically at **http://localhost:8501**.

> If the browser doesn't open, open it manually and go to http://localhost:8501

### Stop the app

Press `Ctrl+C` in the terminal.

### Using the Web App

**Single Mode** — design probes for one gene:
1. Click **Upload FASTA** and select your `.fa` / `.fasta` file, or paste the sequence directly
2. Adjust parameters in the left sidebar (number of probes, oligo length, masking options)
3. Click **Design Probes**
4. View results: probe table, sequence alignment viewer
5. Click download buttons to save `_oligos.txt` and `_seq.txt`

**Batch Mode** — process multiple genes at once:
1. Switch to **Batch** mode at the top of the sidebar
2. Enter the path to a folder containing FASTA files
3. Set an output directory
4. Click **Run Batch** — progress shows per-file
5. Download the batch summary TSV when finished

**Sidebar Parameters:**

| Parameter | Default | What it controls |
|-----------|---------|-----------------|
| Number of probes | 48 | How many probes to design |
| Oligo length | 20 bp | Length of each probe sequence |
| Spacer length | 2 bp | Minimum gap between adjacent probes |
| Target Gibbs FE | −23.0 kcal/mol | Optimal binding free energy |
| Allowable Gibbs range | −26 to −20 | Probes outside this window are excluded |
| Species | human | Which genome/pseudogene databases to use |
| Pseudogene mask | on | Filter probes matching known pseudogenes |
| Genome mask | on | Filter probes in repetitive genome regions |
| RepeatMasker | off | Auto-detect repeats (needs extra installation) |

---

## 6. Using the Command-Line Interface

If you prefer the terminal over the web app, the `probedesign` command is available after activating the environment.

### Top-level commands

```
Usage: probedesign [OPTIONS] COMMAND [ARGS]...

  ProbeDesign - Design oligonucleotide probes for FISH experiments.

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  analyze  Analyze a sequence for probe design feasibility.
  design   Design probes for a target sequence.
```

### `probedesign design` — all options

```
Usage: probedesign design [OPTIONS] INPUT_FILE

  Design probes for a target sequence.

  INPUT_FILE is a FASTA file containing the target sequence.

Options:
  -n, --probes INTEGER            Number of probes to design (default: 48)
  -l, --oligo-length INTEGER      Length of each oligonucleotide (default: 20)
  -s, --spacer-length INTEGER     Minimum gap between probes (default: 2)
  -g, --target-gibbs FLOAT        Target Gibbs free energy in kcal/mol
                                  (default: -23)
  --allowable-gibbs TEXT          Allowable Gibbs FE range as min,max
                                  (default: -26,-20)
  -o, --output TEXT               Output file prefix (default: derived from
                                  input filename)
  -q, --quiet                     Suppress output to stdout
  --species [human|mouse|elegans|drosophila|rat]
                                  Species for masking databases (default:
                                  human)
  --pseudogene-mask / --no-pseudogene-mask
                                  Mask regions that align to pseudogenes
                                  (default: off)
  --genome-mask / --no-genome-mask
                                  Mask repetitive genomic regions (default:
                                  off)
  --index-dir DIRECTORY           Directory containing bowtie indexes
                                  (default: auto-detect)
  --repeatmask-file PATH          FASTA file with N's marking repeat regions
                                  (for manual repeat masking)
  --repeatmask / --no-repeatmask  Run RepeatMasker to automatically mask
                                  repeat regions (default: off)
  --save-bowtie-raw               Save raw bowtie alignment output (large
                                  files)
  --help                          Show this message and exit.
```

### `probedesign analyze` — check a sequence before designing

```
Usage: probedesign analyze [OPTIONS] INPUT_FILE

  Analyze a sequence for probe design feasibility.
  Shows the distribution of Gibbs free energies across the sequence.

Options:
  -l, --oligo-length INTEGER  Length of oligonucleotide to analyze (default: 20)
  --help                      Show this message and exit.
```

### CLI examples

```bash
# Minimal: 32 probes, no masking
probedesign design input.fa --probes 32

# With pseudogene + genome masking (recommended)
probedesign design input.fa --probes 32 --pseudogene-mask --genome-mask

# With a pre-masked FASTA file (N's mark repeat regions)
probedesign design input.fa --probes 32 --repeatmask-file input_masked.fa

# HCR probes (52 bp oligos, different thermodynamic targets)
probedesign design input.fa --probes 20 -l 52 \
  --target-gibbs -60 --allowable-gibbs -80,-40

# Mouse gene, all masking, custom output name
probedesign design my_gene.fa --probes 48 --species mouse \
  --pseudogene-mask --genome-mask -o MyGene_probes

# Check Gibbs FE distribution before designing
probedesign analyze input.fa
```

---

## 7. Masking Options

Masking removes probe positions that would hybridize non-specifically. All masking is optional, but using pseudogene + genome masking is strongly recommended for in vivo experiments.

### Repeat Masking (R Mask) — filters repetitive elements

Three ways to use it, in order of convenience:

| Method | Flag | When to use |
|--------|------|-------------|
| Auto-detect | *(none)* | Your FASTA already has `N`s at repeat positions |
| Manual file | `--repeatmask-file PATH` | You have a pre-masked version of your sequence |
| Automatic | `--repeatmask` | You want RepeatMasker to run automatically (requires extra installation — see [REPEATMASKER.md](../REPEATMASKER.md)) |

`--repeatmask` and `--repeatmask-file` cannot be used together.

### Pseudogene Masking (P Mask) — filters pseudogene cross-hybridization

```bash
probedesign design input.fa --pseudogene-mask
```

Uses 16-mer exact alignment against a curated pseudogene database. Indices are built locally by `setup_all.sh`.

### Genome Masking (B Mask) — filters repetitive genomic regions

```bash
probedesign design input.fa --genome-mask
```

Runs three bowtie searches (12/14/16-mers) against the whole genome. Any position appearing too many times is masked. Genome indices are downloaded by `setup_all.sh`.

### Supported Species

| Species | `--species` | Pseudogene index built? | Genome index downloaded? |
|---------|-------------|------------------------|-------------------------|
| Human | `human` | ✅ by setup_all.sh | ✅ by setup_all.sh |
| Mouse | `mouse` | ✅ by setup_all.sh | ✅ by setup_all.sh |
| Drosophila | `drosophila` | ✅ by setup_all.sh | ✅ by setup_all.sh |
| C. elegans | `elegans` | Manual build needed | Manual build needed |
| Rat | `rat` | Manual build needed | Manual build needed |

---

## 8. Output Files

### `<name>_oligos.txt`

Tab-separated file, one probe per line. Order the sequences in the last two columns from this file.

| Column | Content | Example |
|--------|---------|---------|
| 1 | Probe index (1-based) | `1` |
| 2 | GC content (%) | `55` |
| 3 | Melting temperature Tm (°C) | `66.2` |
| 4 | Gibbs free energy ΔG (kcal/mol) | `−24.7` |
| 5 | Probe sequence (5'→3') | `gtctcagaagctgcgattcg` |
| 6 | Probe name | `MyGene_1` |

### `<name>_seq.txt`

A visual map of the target sequence showing where probes are placed and what was masked. Useful for checking masking behavior.

The file is wrapped at 110 characters per line, with each block showing:
- **Line 1**: Original sequence
- **R line** (if repeat masking used): `R` = masked, base = unmasked
- **P line** (if pseudogene masking used): `P` = masked
- **B line** (if genome masking used): `B` = masked
- **F line** (always): `F` = thermodynamically excluded, base = valid
- **Probe lines**: complementary sequence + probe number, ΔG, GC%

### Additional files (masking enabled)

| File | Content |
|------|---------|
| `<name>_bowtie_pseudogene.txt` | Per-position hit counts from pseudogene search |
| `<name>_bowtie_genome.txt` | Per-position hit counts from genome search |
| `<name>_bowtie_*_raw.txt` | Raw bowtie alignment output (only with `--save-bowtie-raw`) |

---

## 9. Manual Environment Setup (Advanced)

If you prefer to create the environment yourself instead of using `setup_all.sh`:

### Using environment.yml (recommended manual method)

```bash
# From the repository root:

# micromamba:
micromamba env create -f streamlit_app/environment.yml

# mamba:
mamba env create -f streamlit_app/environment.yml

# conda:
conda env create -f streamlit_app/environment.yml

# Then install the probedesign package:
micromamba activate probedesign   # (or mamba / conda)
pip install -e .
```

### Without environment.yml

```bash
micromamba create -n probedesign \
    -c conda-forge -c bioconda \
    --channel-priority strict \
    python=3.11 "click>=8.0" bowtie=1.3.1 "pytest>=7.0" -y

micromamba activate probedesign
pip install -e .
pip install "streamlit>=1.32" "pandas>=1.5"
```

### Updating an existing environment

```bash
micromamba env update -n probedesign \
    -f streamlit_app/environment.yml --prune
```

### Removing the environment (clean slate)

```bash
micromamba env remove -n probedesign
# (or: mamba / conda env remove -n probedesign)
```

---

## 10. Troubleshooting

### "micromamba / mamba / conda: command not found"

The package manager is not on your PATH. Close and reopen the terminal (the installer modifies your shell profile). If that doesn't help, re-run the package manager installer.

### "bowtie: command not found" after activating the environment

Make sure you are using the probedesign environment:
```bash
micromamba activate probedesign
bowtie --version   # should print: bowtie-align-s version 1.3.1
```

If bowtie shows a **version mismatch** or is the wrong tool (e.g., a file manager), you may have the wrong bowtie installed. On macOS, `brew install bowtie` installs an unrelated app. Only install bowtie via **bioconda**:
```bash
micromamba install -n probedesign -c bioconda bowtie=1.3.1
```

### Setup script fails at "Build pseudogene indices"

Check that the pseudogene FASTA files exist:
```bash
ls probedesign/pseudogeneDBs/*.fasta
# Should list: human.fasta  mouse.fasta  drosophila.fasta
```
If missing, re-clone the repository.

### Genome download interrupted (partial download)

The script re-downloads if fewer than 6 `.bt2` files are present for a species. Simply re-run:
```bash
./streamlit_app/setup_all.sh
```

### "streamlit: command not found"

Streamlit is installed into the `probedesign` environment. Make sure it is activated:
```bash
micromamba activate probedesign
streamlit run streamlit_app/app.py
```

### Port 8501 already in use

Another Streamlit instance is running. Either stop it (Ctrl+C in the other terminal) or use a different port:
```bash
streamlit run streamlit_app/app.py --server.port 8502
```

### Windows: script fails with "bash: ./streamlit_app/setup_all.sh: No such file or directory"

Make sure you are inside the WSL2 Ubuntu terminal (not PowerShell or CMD). Navigate to the repository inside WSL2:
```bash
# Example — adjust the path to match where you cloned the repo on Windows:
cd /mnt/c/Users/YourName/Documents/ProbeDesign
./streamlit_app/setup_all.sh
```

---

## 11. Further Reading

| Document | Description |
|----------|-------------|
| [probe_design_principles.md](probe_design_principles.md) | Complete reference: thermodynamic model, masking algorithm, DP optimization, all parameters |
| [summary.md](summary.md) | Development notes: installation findings, Bowtie 1/2 index compatibility, source code changes |
| [environment.yml](environment.yml) | Conda environment specification (all dependencies with pinned versions) |
| [setup_all.sh](setup_all.sh) | Automated setup script source |
| [../BOWTIE.md](../BOWTIE.md) | Bowtie index building guide |
| [../REPEATMASKER.md](../REPEATMASKER.md) | RepeatMasker installation guide |
