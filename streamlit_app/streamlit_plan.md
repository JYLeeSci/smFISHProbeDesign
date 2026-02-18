# Plan: Streamlit GUI for ProbeDesign
TL;DR: Build a streamlit_app/ directory with app.py as the main entry point. The app wraps design_probes() from core.py directly (no subprocess), with a sidebar for all parameters and a main panel for input/results. Two top-level modes — single and batch — are toggled via a radio button. Results show an oligos table, scrollable _seq.txt viewer, and download buttons. All external-tool checks (bowtie, index files) surface as user-friendly warnings before running.

## Steps

1. Create directory and file scaffold

Create streamlit_app/app.py, streamlit_app/utils.py, streamlit_app/requirements.txt
Add streamlit>=1.32 to requirements.txt (keep it separate from pyproject.toml to avoid adding a heavy dependency to the core library)

2. Build utils.py — backend helpers

fasta_from_paste(text) -> Path: validates pasted FASTA text (checks >header line present, valid bases), writes to a tempfile.NamedTemporaryFile, returns path
run_design(input_path, output_dir, params) -> ProbeDesignResult: sets CWD to output_dir, calls design_probes() from core.py, restores CWD, returns result
check_bowtie() -> bool: runs shutil.which("bowtie") and verifies version string
check_index_dir(index_dir, species, pseudogene, genome) -> List[str]: returns list of missing index file warnings

3. Build app.py — sidebar (parameters panel)

Mode toggle at top: st.radio("Mode", ["Single sequence", "Batch"]) — controls which main panel is shown
Basic parameters section:
    st.number_input for Number of probes (default 48)
    st.number_input for Oligo length (default 20)
    st.number_input for Spacer length (default 2)
    st.number_input for Target Gibbs (default −23.0)
    Two st.number_input for Allowable Gibbs min/max (defaults −26, −20)
    st.selectbox for Species: human, mouse, elegans, drosophila, rat
Masking section (collapsible st.expander):
    st.checkbox for Pseudogene mask
    st.checkbox for Genome mask
    st.radio for Repeat masking: None / Auto (RepeatMasker) / Provide file — enforces mutual exclusion
    Conditional st.file_uploader shown only when "Provide file" is selected
    st.text_input for Bowtie index directory (defaults to workspace bowtie_indexes)
    st.checkbox for Save raw bowtie output (advanced)
Output section:
    st.text_input for Output directory (defaults to ~/Desktop/probedesign_output)
    st.text_input for Output name prefix (optional, auto-filled from filename)

4. Build app.py — Single mode main panel

st.tabs(["Upload FASTA", "Paste FASTA"]) for input method
    Upload tab: st.file_uploader(type=["fa","fasta","txt"])
    Paste tab: st.text_area for full FASTA format text
Prerequisite checks run on st.button("Design Probes") click, before the actual run:
    If masking enabled: check bowtie present and index files exist → st.warning per issue
st.spinner("Running probe design...") wraps the run_design() call
After success: display results (step 6)
After failure: st.error with traceback summary

5. Build app.py — Batch mode main panel

st.text_input for Input directory path (must contain .fa/.fasta files)
st.button("Run Batch") — iterates FASTA files one-by-one:
    st.progress bar tracks overall progress
    Per-file status shown in st.status (expandable) showing ✓ / ✗ per file
After all runs: summary table of file → probe count → score → status
Download all outputs as a .zip via st.download_button

6. Results display (shared between single and batch)

Oligos table: st.dataframe with columns Index, GC%, Tm, Gibbs, Sequence, Name — with column sorting enabled
Sequence visualization: st.text (monospaced) inside st.expander("Sequence alignment (_seq.txt)") — scrollable pre-formatted text
Download buttons:
    st.download_button for _oligos.txt
    st.download_button for _seq.txt
    st.download_button for parsed bowtie_genome and/or pseudogene hits .txt
    If bowtie raw saved: additional buttons for raw files

7. Output directory handling

utils.run_design() creates output_dir if it doesn't exist (os.makedirs)
Uses os.chdir(output_dir) + restore pattern so write_output_files() in output.py writes to the correct location

8. Add launch instructions — streamlit_app/README.md:

    micromamba activate probedesign
    pip install streamlit
    streamlit run streamlit_app/app.py

