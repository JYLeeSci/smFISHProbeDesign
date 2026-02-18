# ProbeDesign Streamlit App

Web-based GUI for designing oligonucleotide probes for single molecule RNA FISH.

## Installation

The app requires streamlit, which should be installed in your probedesign environment:

```bash
micromamba activate probedesign
pip install streamlit
```

## Launch

From the repository root:

```bash
micromamba activate probedesign
streamlit run streamlit_app/app.py
```

Or from within the `streamlit_app/` directory:

```bash
micromamba activate probedesign
streamlit run app.py
```

The app will open in your browser at http://localhost:8501

## Features

### Single Mode
- Upload FASTA files or paste sequences directly
- Configure all probe design parameters via sidebar
- View designed probes in an interactive table
- Download oligos and sequence alignment files
- Support for pseudogene, genome, and repeat masking

### Batch Mode
- Process multiple FASTA files at once
- Input from local directory or file upload
- Progress tracking during batch processing
- Batch summary with probe counts and scores per file
- All outputs written to specified directory

## Requirements

- Python 3.8+
- probedesign package (installed in the same environment)
- streamlit>=1.32

Optional (for masking):
- bowtie (via conda: `mamba install bowtie`)
- RepeatMasker (via conda: `mamba install -c bioconda -c conda-forge repeatmasker`)

See [BOWTIE.md](../BOWTIE.md) and [REPEATMASKER.md](../REPEATMASKER.md) in the repository root for setup instructions.
