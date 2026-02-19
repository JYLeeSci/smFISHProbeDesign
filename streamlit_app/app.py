"""Streamlit GUI for ProbeDesign — single molecule RNA FISH probe designer."""

import os
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

# Ensure the probedesign package is importable from the repo root
_repo_root = Path(__file__).resolve().parent.parent
if str(_repo_root / "src") not in sys.path:
    sys.path.insert(0, str(_repo_root / "src"))

from utils import (
    validate_fasta_text,
    save_fasta_to_temp,
    save_uploaded_file_to_temp,
    cleanup_temp_files,
    check_prerequisites,
    run_design,
    run_batch,
    discover_fasta_files,
    format_batch_summary,
    write_batch_summary,
    package_batch_results_zip,
    BatchResult,
    clean_name_for_prefix,
    extract_first_fasta_header,
)
from probedesign.output import (
    format_oligos_content,
    format_seq_content,
    format_hits_content,
)


# ── Page config ──────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="smFISH ProbeDesign",
    layout="wide",
)

# ── Session state initialization ─────────────────────────────────────────────

for key, default in {
    "single_result": None,
    "single_input_name": None,
    "batch_results": None,
    "batch_output_dir": None,
    "temp_files": [],
}.items():
    if key not in st.session_state:
        st.session_state[key] = default


# ── Sidebar — parameters panel ───────────────────────────────────────────────

with st.sidebar:
    st.header("Probe Design")

    mode = st.radio(
        "Mode",
        ["Single sequence", "Batch"],
        key="mode_radio",
        horizontal=True,
    )

    st.write("")

    # Basic parameters
    n_probes = st.number_input(
        "Number of probes",
        min_value=1, max_value=200, value=48, step=1,
        key="n_probes",
    )
    oligo_mode = st.radio(
        "Oligo length mode",
        ["Fixed", "Mixed range"],
        key="oligo_mode",
        horizontal=True,
        help="Fixed: all probes same length. Mixed range: probes vary between min-max length.",
    )

    mixed_lengths = None
    if oligo_mode == "Fixed":
        oligo_length = st.number_input(
            "Oligo length (bp)",
            min_value=10, max_value=60, value=20, step=1,
            key="oligo_length",
        )
    else:
        ocol1, ocol2 = st.columns(2)
        with ocol1:
            mixed_min = st.number_input(
                "Min length (bp)",
                min_value=10, max_value=60, value=18, step=1,
                key="mixed_min",
            )
        with ocol2:
            mixed_max = st.number_input(
                "Max length (bp)",
                min_value=10, max_value=60, value=22, step=1,
                key="mixed_max",
            )
        if mixed_min >= mixed_max:
            st.warning("Min length must be less than max length.")
        oligo_length = mixed_min
        mixed_lengths = (mixed_min, mixed_max)
    spacer_length = st.number_input(
        "Spacer length (bp)",
        min_value=0, max_value=20, value=2, step=1,
        key="spacer_length",
    )
    target_gibbs = st.number_input(
        "Target Gibbs FE (kcal/mol)",
        min_value=-100.0, max_value=0.0, value=-23.0, step=0.5,
        format="%.1f",
        key="target_gibbs",
    )

    col1, col2 = st.columns(2)
    with col1:
        gibbs_min = st.number_input(
            "Gibbs min",
            min_value=-200.0, max_value=0.0, value=-26.0, step=0.5,
            format="%.1f",
            key="gibbs_min",
        )
    with col2:
        gibbs_max = st.number_input(
            "Gibbs max",
            min_value=-200.0, max_value=0.0, value=-20.0, step=0.5,
            format="%.1f",
            key="gibbs_max",
        )

    species = st.selectbox(
        "Species",
        ["human", "mouse", "elegans", "drosophila", "rat"],
        key="species",
    )

    st.divider()

    # Bowtie masking — prominent, on by default
    pseudogene_mask_on = st.checkbox(
        "Pseudogene mask", value=True, key="pseudogene_mask",
        help="Masking against pseudogenes"
    )
    genome_mask_on = st.checkbox(
        "Genome mask", value=True, key="genome_mask",
        help="Masking against the whole genome (more comprehensive but slower than pseudogene masking"
    )

    index_dir = None
    if pseudogene_mask_on or genome_mask_on:
        index_dir = st.text_input(
            "Bowtie index directory",
            value=str(_repo_root / "bowtie_indexes"),
            key="index_dir",
            help="Directory containing bowtie indexes for the selected species (required if masking is enabled)",
        )

    save_bowtie_raw = st.checkbox(
        "Save raw bowtie output",
        value=False,
        key="save_bowtie_raw",
        help="Include full bowtie alignment files in downloads",
    )

    st.divider()

    # RepeatMasker — advanced, hidden by default
    with st.expander("Advanced RepeatMasker options", expanded=False):
        repeatmask_mode = st.radio(
            "Repeat masking",
            ["None", "Auto (RepeatMasker)", "Provide file"],
            key="repeatmask_mode",
        )

        repeatmask_upload = None
        if repeatmask_mode == "Provide file":
            repeatmask_upload = st.file_uploader(
                "Repeat mask FASTA",
                type=["fa", "fasta", "txt"],
                key="repeatmask_file_upload",
            )


# ── Helper: map repeatmask_mode UI string to internal key ────────────────────

_RM_MODE_MAP = {
    "None": "none",
    "Auto (RepeatMasker)": "auto",
    "Provide file": "file",
}
rm_mode_key = _RM_MODE_MAP[repeatmask_mode]


# ── Helper: run prerequisite checks and display warnings ─────────────────────

def _run_prereq_checks() -> bool:
    """Run prerequisite checks. Returns True if OK to proceed."""
    prereq = check_prerequisites(
        pseudogene_mask=pseudogene_mask_on,
        genome_mask=genome_mask_on,
        repeatmask_mode=rm_mode_key,
        species=species,
        index_dir=index_dir or "",
    )

    ok = True

    if (pseudogene_mask_on or genome_mask_on) and not prereq.bowtie_ok:
        st.error("Bowtie is not installed or not found. Masking requires bowtie.")
        ok = False

    if rm_mode_key == "auto" and not prereq.repeatmasker_ok:
        st.error("RepeatMasker is not installed. Auto repeat masking requires RepeatMasker.")
        ok = False

    for msg in prereq.missing_indexes:
        st.warning(msg)

    return ok


# ── Helper: resolve repeat mask file from upload ─────────────────────────────

def _resolve_repeatmask_file():
    """Resolve the repeat mask file path. Returns None if not applicable."""
    if rm_mode_key == "file":
        if repeatmask_upload is None:
            st.warning("Please upload a repeat mask file.")
            return "__STOP__"
        tmp = save_uploaded_file_to_temp(repeatmask_upload)
        st.session_state["temp_files"].append(tmp)
        return str(tmp)
    return None


# ── Helper: display single-mode results ──────────────────────────────────────

def _show_single_results():
    """Display results stored in session state for single mode."""
    run_result = st.session_state.get("single_result")
    if run_result is None:
        return

    if run_result.error:
        st.error(run_result.error)

    # Show design log
    if run_result.stdout_log.strip():
        with st.expander("Design log"):
            st.code(run_result.stdout_log, language=None)

    if run_result.result is not None and run_result.result.probes:
        result = run_result.result
        input_name = st.session_state.get("single_input_name") or "probes"
        prefix = Path(input_name).stem

        # Summary metrics
        mcol1, mcol2 = st.columns(2)
        mcol1.metric("Probes found", len(result.probes))
        mcol2.metric("Score", f"{result.score:.4f}")

        # Probes table
        st.subheader("Designed Probes")
        probe_data = []
        probe_lengths = set(len(p.sequence) for p in result.probes)
        show_length = len(probe_lengths) > 1
        for p in result.probes:
            row = {"Index": p.index}
            if show_length:
                row["Length"] = len(p.sequence)
            row.update({
                "GC%": p.gc_percent,
                "Tm": p.tm,
                "Gibbs FE": p.gibbs_fe,
                "Sequence": p.sequence,
                "Name": p.name,
            })
            probe_data.append(row)
        df = pd.DataFrame(probe_data)
        st.dataframe(df, width='stretch', hide_index=True)

        # Sequence alignment viewer
        seq_content = format_seq_content(result, mask_seqs=result.mask_strings)
        with st.expander("Sequence alignment"):
            st.code(seq_content, language=None)

        # Download buttons
        st.subheader("Download Files")
        oligos_content = format_oligos_content(result)

        dcol1, dcol2 = st.columns(2)
        with dcol1:
            st.download_button(
                label="Download oligos.txt",
                data=oligos_content,
                file_name=f"{prefix}_oligos.txt",
                mime="text/plain",
                key="dl_oligos",
            )
        with dcol2:
            st.download_button(
                label="Download seq.txt",
                data=seq_content,
                file_name=f"{prefix}_seq.txt",
                mime="text/plain",
                key="dl_seq",
            )

        # Bowtie hits downloads
        hits_files = format_hits_content(result)
        for suffix, content in hits_files.items():
            st.download_button(
                label=f"Download {suffix}",
                data=content,
                file_name=f"{prefix}_{suffix}",
                mime="text/plain",
                key=f"dl_{suffix}",
            )

        # Raw bowtie downloads
        if result.bowtie_pseudogene_raw:
            st.download_button(
                label="Download bowtie_pseudogene_raw.txt",
                data=result.bowtie_pseudogene_raw,
                file_name=f"{prefix}_bowtie_pseudogene_raw.txt",
                mime="text/plain",
                key="dl_pseudo_raw",
            )
        if result.bowtie_genome_raw:
            st.download_button(
                label="Download bowtie_genome_raw.txt",
                data=result.bowtie_genome_raw,
                file_name=f"{prefix}_bowtie_genome_raw.txt",
                mime="text/plain",
                key="dl_genome_raw",
            )


# ══════════════════════════════════════════════════════════════════════════════
#  SINGLE MODE
# ══════════════════════════════════════════════════════════════════════════════

if mode == "Single sequence":
    st.title("ProbeDesign")
    st.caption("Design oligonucleotide probes for single molecule RNA FISH")

    # Output name (shown above tabs so it's always visible)
    output_name = st.text_input(
        "Output name prefix (optional, auto-derived from filename)",
        value="",
        key="single_output_name",
    )

    # Gibbs range validation
    if gibbs_min >= gibbs_max:
        st.warning("Gibbs min should be less than Gibbs max.")

    # Input tabs — each tab has its own run button to unambiguously signal input source
    upload_tab, paste_tab = st.tabs(["Upload FASTA", "Paste FASTA"])

    with upload_tab:
        uploaded_file = st.file_uploader(
            "Choose a FASTA file",
            type=["fa", "fasta", "txt"],
            key="single_upload",
        )
        run_from_upload = st.button("Design Probes", key="single_run_upload", type="primary")

    with paste_tab:
        pasted_text = st.text_area(
            "Paste FASTA-formatted sequence",
            height=200,
            placeholder=">gene_name\nATCGATCG...",
            key="single_paste",
        )
        run_from_paste = st.button("Design Probes", key="single_run_paste", type="primary")

    if run_from_upload or run_from_paste:
        # Clean up previous temp files
        cleanup_temp_files(st.session_state["temp_files"])
        st.session_state["temp_files"] = []
        st.session_state["single_result"] = None

        input_path = None
        input_name = None
        auto_name = None

        if run_from_upload:
            if uploaded_file is None:
                st.warning("Please upload a FASTA file.")
                st.stop()
            tmp_path = save_uploaded_file_to_temp(uploaded_file)
            st.session_state["temp_files"].append(tmp_path)
            input_path = str(tmp_path)
            input_name = uploaded_file.name
            auto_name = clean_name_for_prefix(Path(uploaded_file.name).stem)
        else:  # run_from_paste
            if not pasted_text.strip():
                st.warning("Please paste a FASTA sequence.")
                st.stop()
            valid, err = validate_fasta_text(pasted_text)
            if not valid:
                st.error(f"Invalid FASTA: {err}")
                st.stop()
            tmp_path = save_fasta_to_temp(pasted_text)
            st.session_state["temp_files"].append(tmp_path)
            input_path = str(tmp_path)
            input_name = "pasted_sequence"
            header = extract_first_fasta_header(pasted_text)
            auto_name = clean_name_for_prefix(header) if header else "probe"

        # Resolve repeatmask file
        rm_file = _resolve_repeatmask_file()
        if rm_file == "__STOP__":
            st.stop()

        # Prerequisite checks
        if not _run_prereq_checks():
            st.stop()

        # Run design — prefer explicit user input, fall back to auto-derived name
        resolved_name = output_name.strip() if output_name.strip() else auto_name

        with st.spinner("Running probe design..."):
            run_result = run_design(
                input_path=input_path,
                n_probes=n_probes,
                oligo_length=oligo_length,
                spacer_length=spacer_length,
                target_gibbs=target_gibbs,
                allowable_gibbs=(gibbs_min, gibbs_max),
                output_name=resolved_name,
                species=species,
                pseudogene_mask=pseudogene_mask_on,
                genome_mask=genome_mask_on,
                index_dir=index_dir,
                repeatmask_mode=rm_mode_key,
                repeatmask_file=rm_file,
                save_bowtie_raw=save_bowtie_raw,
                mixed_lengths=mixed_lengths,
            )

        st.session_state["single_result"] = run_result
        st.session_state["single_input_name"] = input_name

    # Always render results if available
    _show_single_results()


# ══════════════════════════════════════════════════════════════════════════════
#  BATCH MODE
# ══════════════════════════════════════════════════════════════════════════════

else:
    st.title("ProbeDesign — Batch Mode")
    st.caption("Design probes for multiple sequences at once")

    # Input method
    input_method = st.radio(
        "Input method",
        ["Upload files", "Directory path"],
        key="batch_input_method",
        horizontal=True,
    )

    input_dir = None
    uploaded_files = None

    if input_method == "Directory path":
        input_dir = st.text_input(
            "Directory containing .fa/.fasta files",
            key="batch_input_dir",
        )
        st.caption("Input filesystem path. Use 'Upload files' if connecting remotely.")
    else:
        uploaded_files = st.file_uploader(
            "Upload FASTA files",
            type=["fa", "fasta", "txt"],
            accept_multiple_files=True,
            key="batch_upload",
        )

    # Output directory — optional, for server-side power users
    with st.expander("Auto-save outputs to Local/Server directory (optional)", expanded=False):
        st.caption("Leave blank to download results as a ZIP file (recommended for remote access).")
        batch_output_dir = st.text_input(
            "Output directory path on server",
            value="",
            key="batch_output_dir_input",
        )
    
    st.write("")

    # Gibbs range validation
    if gibbs_min >= gibbs_max:
        st.warning("Gibbs min should be less than Gibbs max.")

    # Run button
    if st.button("Run Batch", key="batch_run_btn", type="primary"):
        # Clean up previous temp files
        cleanup_temp_files(st.session_state["temp_files"])
        st.session_state["temp_files"] = []
        st.session_state["batch_results"] = None
        st.session_state["batch_output_dir"] = None

        # Resolve input files
        fasta_paths = []
        batch_name_overrides = {}
        if input_method == "Directory path":
            if not input_dir or not Path(input_dir).is_dir():
                st.error(f"Directory not found: {input_dir}")
                st.stop()
            fasta_paths = discover_fasta_files(input_dir)
            if not fasta_paths:
                st.warning("No .fa or .fasta files found in directory.")
                st.stop()
        else:
            if not uploaded_files:
                st.warning("Please upload at least one FASTA file.")
                st.stop()
            for uf in uploaded_files:
                tmp = save_uploaded_file_to_temp(uf)
                st.session_state["temp_files"].append(tmp)
                fasta_paths.append(tmp)
                batch_name_overrides[tmp] = clean_name_for_prefix(Path(uf.name).stem)

        # Resolve repeatmask file
        rm_file = _resolve_repeatmask_file()
        if rm_file == "__STOP__":
            st.stop()

        # Prerequisite checks
        if not _run_prereq_checks():
            st.stop()

        # Build params dict
        params = {
            "n_probes": n_probes,
            "oligo_length": oligo_length,
            "spacer_length": spacer_length,
            "target_gibbs": target_gibbs,
            "allowable_gibbs": (gibbs_min, gibbs_max),
            "species": species,
            "pseudogene_mask": pseudogene_mask_on,
            "genome_mask": genome_mask_on,
            "index_dir": index_dir,
            "repeatmask_mode": rm_mode_key,
            "repeatmask_file": rm_file,
            "save_bowtie_raw": save_bowtie_raw,
            "mixed_lengths": mixed_lengths,
        }

        # Run batch with progress
        progress_bar = st.progress(0, text="Starting batch...")

        def _progress_cb(current, total, filename):
            if total > 0:
                progress_bar.progress(
                    current / total,
                    text=f"Processing {filename} ({current + 1}/{total})..."
                    if current < total
                    else "Batch complete!",
                )

        batch_results = run_batch(
            fasta_paths=fasta_paths,
            output_dir=batch_output_dir.strip() or None,
            params=params,
            progress_callback=_progress_cb,
            name_overrides=batch_name_overrides if batch_name_overrides else None,
        )

        # Write summary file to disk only if a server directory was specified
        if batch_output_dir.strip():
            write_batch_summary(batch_results, batch_output_dir)

        progress_bar.progress(1.0, text="Batch complete!")
        st.session_state["batch_results"] = batch_results
        st.session_state["batch_output_dir"] = batch_output_dir.strip() or None

    # Always render batch results if available
    batch_results = st.session_state.get("batch_results")
    saved_output_dir = st.session_state.get("batch_output_dir")

    if batch_results is not None:
        st.subheader("Batch Summary")

        # Summary table
        summary_data = []
        for br in batch_results:
            summary_data.append({
                "File": br.filename,
                "Probes": br.n_probes_found,
                "Score": f"{br.score:.4f}" if br.score != float('inf') else "N/A",
                "Status": br.status,
                "Error": br.error_msg or "",
            })
        df = pd.DataFrame(summary_data)
        st.dataframe(df, width='stretch', hide_index=True)

        # Summary download buttons
        zip_bytes = package_batch_results_zip(batch_results)
        dcol1, dcol2 = st.columns(2)
        with dcol1:
            st.download_button(
                "Download all results (.zip)",
                data=zip_bytes,
                file_name="probedesign_batch_results.zip",
                mime="application/zip",
                key="dl_batch_zip",
            )
        with dcol2:
            summary_text = format_batch_summary(batch_results)
            st.download_button(
                "Download batch_summary.tsv",
                data=summary_text,
                file_name="batch_summary.tsv",
                mime="text/tab-separated-values",
                key="dl_batch_summary",
            )

        if saved_output_dir:
            st.info(f"Output files also written to server: {saved_output_dir}")

        # Per-file expandable details
        for idx, br in enumerate(batch_results):
            if br.result is not None and br.result.probes:
                with st.expander(f"{br.filename} — {br.n_probes_found} probes (score: {br.score:.4f})"):
                    br_lengths = set(len(p.sequence) for p in br.result.probes)
                    br_show_length = len(br_lengths) > 1
                    probe_data = []
                    for p in br.result.probes:
                        row = {"Index": p.index}
                        if br_show_length:
                            row["Length"] = len(p.sequence)
                        row.update({
                            "GC%": p.gc_percent,
                            "Tm": p.tm,
                            "Gibbs FE": p.gibbs_fe,
                            "Sequence": p.sequence,
                        })
                        probe_data.append(row)
                    st.dataframe(
                        pd.DataFrame(probe_data),
                        width='stretch',
                        hide_index=True,
                    )
            elif br.error_msg:
                with st.expander(f"{br.filename} — FAILED"):
                    st.error(br.error_msg)
