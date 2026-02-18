"""Backend helpers for the Streamlit ProbeDesign app."""

import io
import os
import re
import tempfile
from contextlib import redirect_stdout
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from probedesign.core import design_probes, ProbeDesignResult
from probedesign.output import (
    write_output_files,
    format_oligos_content,
    format_seq_content,
    format_hits_content,
)
from probedesign.masking import find_bowtie, find_repeatmasker, run_repeatmasker


# --- Species-to-index maps (mirrors masking.py) ---

PSEUDO_DB = {
    "human": "humanPseudo",
    "mouse": "mousePseudo",
    "elegans": "celegansPseudo",
    "drosophila": "drosophilaPseudo",
    "rat": "ratPseudo",
}

GENOME_DB = {
    "human": "GCA_000001405.15_GRCh38_no_alt_analysis_set",
    "mouse": "mm10",
    "elegans": "celegans",
    "drosophila": "drosophila",
    "rat": "rat",
}


# --- FASTA validation and temp file helpers ---


def validate_fasta_text(text: str) -> Tuple[bool, str]:
    """Validate pasted FASTA text.

    Checks:
    - At least one header line starting with '>'
    - At least one non-empty sequence line after each header
    - Sequence lines contain only valid bases (ACGTNacgtn and whitespace)

    Returns:
        (is_valid, error_message). error_message is "" on success.
    """
    text = text.strip()
    if not text:
        return False, "Empty input."

    lines = text.split("\n")

    # Must start with a header
    if not lines[0].strip().startswith(">"):
        return False, "FASTA must start with a header line beginning with '>'."

    found_header = False
    found_seq_after_header = False

    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if found_header and not found_seq_after_header:
                return False, "Header line has no sequence following it."
            found_header = True
            found_seq_after_header = False
        else:
            if not found_header:
                return False, "Sequence data found before any header line."
            # Check for valid characters (ACGTN only, case-insensitive)
            if not re.match(r'^[ACGTNacgtn]+$', line):
                invalid = re.findall(r'[^ACGTNacgtn]', line)
                return False, f"Invalid characters in sequence: {set(invalid)}"
            found_seq_after_header = True

    if not found_seq_after_header:
        return False, "No sequence data found after the last header."

    return True, ""


def save_fasta_to_temp(text: str) -> Path:
    """Write validated FASTA text to a named temp file.

    The file is NOT auto-deleted — the caller is responsible for cleanup.

    Returns:
        Path to the temp file.
    """
    f = tempfile.NamedTemporaryFile(
        mode='w', suffix='.fa', delete=False, prefix='probedesign_'
    )
    f.write(text)
    f.close()
    return Path(f.name)


def save_uploaded_file_to_temp(uploaded_file) -> Path:
    """Write a Streamlit UploadedFile to a named temp file.

    Returns:
        Path to the temp file.
    """
    content = uploaded_file.getvalue().decode('utf-8')
    f = tempfile.NamedTemporaryFile(
        mode='w', suffix='.fa', delete=False, prefix='probedesign_'
    )
    f.write(content)
    f.close()
    return Path(f.name)


def cleanup_temp_files(paths: List[Path]) -> None:
    """Delete temporary files. Silently ignores missing files."""
    for p in paths:
        try:
            if p.exists():
                p.unlink()
        except OSError:
            pass


# --- Prerequisite checks ---


@dataclass
class PrereqStatus:
    bowtie_ok: bool = False
    bowtie_path: Optional[str] = None
    repeatmasker_ok: bool = False
    repeatmasker_path: Optional[str] = None
    index_dir_exists: bool = False
    missing_indexes: List[str] = field(default_factory=list)


def check_prerequisites(
    pseudogene_mask: bool,
    genome_mask: bool,
    repeatmask_mode: str,
    species: str,
    index_dir: str,
) -> PrereqStatus:
    """Check that all required external tools and indexes are available.

    Args:
        pseudogene_mask: Whether pseudogene masking is enabled
        genome_mask: Whether genome masking is enabled
        repeatmask_mode: One of "none", "auto", "file"
        species: Species for masking databases
        index_dir: Path to bowtie indexes directory

    Returns:
        PrereqStatus with detailed availability info
    """
    status = PrereqStatus()

    # Check bowtie
    if pseudogene_mask or genome_mask:
        try:
            status.bowtie_path = find_bowtie()
            status.bowtie_ok = True
        except (FileNotFoundError, Exception):
            status.bowtie_ok = False

    # Check RepeatMasker
    if repeatmask_mode == "auto":
        try:
            status.repeatmasker_path = find_repeatmasker()
            status.repeatmasker_ok = True
        except (FileNotFoundError, Exception):
            status.repeatmasker_ok = False

    # Check index directory and files
    # Bowtie 1.3+ can read both .ebwt (bowtie1) and .bt2 (bowtie2) indexes
    if pseudogene_mask or genome_mask:
        if index_dir and Path(index_dir).is_dir():
            status.index_dir_exists = True

            if pseudogene_mask:
                db_name = PSEUDO_DB.get(species.lower(), "")
                if db_name:
                    d = Path(index_dir)
                    if not (d / f"{db_name}.1.ebwt").exists() and not (d / f"{db_name}.1.bt2").exists():
                        status.missing_indexes.append(
                            f"Pseudogene index not found: {db_name} in {index_dir}"
                        )

            if genome_mask:
                db_name = GENOME_DB.get(species.lower(), "")
                if db_name:
                    d = Path(index_dir)
                    if not (d / f"{db_name}.1.ebwt").exists() and not (d / f"{db_name}.1.bt2").exists():
                        status.missing_indexes.append(
                            f"Genome index not found: {db_name} in {index_dir}"
                        )
        else:
            status.index_dir_exists = False
            if index_dir:
                status.missing_indexes.append(
                    f"Index directory does not exist: {index_dir}"
                )
            else:
                status.missing_indexes.append(
                    "No index directory specified. Required for pseudogene/genome masking."
                )

    return status


# --- Design runner ---


@dataclass
class DesignRunResult:
    result: Optional[ProbeDesignResult] = None
    stdout_log: str = ""
    error: Optional[str] = None


def run_design(
    input_path: str,
    n_probes: int = 48,
    oligo_length: int = 20,
    spacer_length: int = 2,
    target_gibbs: float = -23.0,
    allowable_gibbs: Tuple[float, float] = (-26.0, -20.0),
    output_name: Optional[str] = None,
    species: str = "human",
    pseudogene_mask: bool = False,
    genome_mask: bool = False,
    index_dir: Optional[str] = None,
    repeatmask_mode: str = "none",
    repeatmask_file: Optional[str] = None,
    save_bowtie_raw: bool = False,
) -> DesignRunResult:
    """Run probe design, replicating the CLI orchestration logic.

    1. If repeatmask_mode == "auto", calls run_repeatmasker() first.
    2. Calls design_probes() with the resolved repeatmask_file.
    3. Captures stdout (design_probes prints progress messages).
    4. Catches all exceptions and returns them as error string.

    Returns:
        DesignRunResult with either result or error populated.
    """
    stdout_capture = io.StringIO()
    try:
        actual_repeatmask_file = repeatmask_file

        # Step 1: RepeatMasker (mirrors cli.py lines 144-162)
        if repeatmask_mode == "auto":
            with redirect_stdout(stdout_capture):
                print(f"Running RepeatMasker on {input_path}...")
            masked_file = run_repeatmasker(Path(input_path), species=species)
            actual_repeatmask_file = str(masked_file)
            with redirect_stdout(stdout_capture):
                print(f"RepeatMasker output: {masked_file}")

        # Step 2: design with stdout capture
        with redirect_stdout(stdout_capture):
            result = design_probes(
                input_file=input_path,
                n_probes=n_probes,
                oligo_length=oligo_length,
                spacer_length=spacer_length,
                target_gibbs=target_gibbs,
                allowable_gibbs=allowable_gibbs,
                output_name=output_name,
                species=species,
                pseudogene_mask=pseudogene_mask,
                genome_mask=genome_mask,
                index_dir=index_dir,
                repeatmask_file=actual_repeatmask_file,
                save_bowtie_raw=save_bowtie_raw,
            )

        if not result.probes:
            return DesignRunResult(
                result=result,
                stdout_log=stdout_capture.getvalue(),
                error="No valid probes found for this sequence. "
                      "Try adjusting Gibbs FE range or oligo length.",
            )

        return DesignRunResult(
            result=result,
            stdout_log=stdout_capture.getvalue(),
            error=None,
        )

    except Exception as e:
        return DesignRunResult(
            result=None,
            stdout_log=stdout_capture.getvalue(),
            error=f"{type(e).__name__}: {e}",
        )


# --- Batch runner ---


@dataclass
class BatchResult:
    filename: str
    n_probes_found: int = 0
    score: float = float('inf')
    status: str = "pending"  # "success", "error", "no_probes"
    error_msg: Optional[str] = None
    result: Optional[ProbeDesignResult] = None


def discover_fasta_files(directory: str) -> List[Path]:
    """Find all .fa/.fasta files in a directory (non-recursive).

    Returns:
        Sorted list of FASTA file paths.
    """
    d = Path(directory)
    files = list(d.glob("*.fa")) + list(d.glob("*.fasta"))
    return sorted(set(files), key=lambda p: p.name)


def run_batch(
    fasta_paths: List[Path],
    output_dir: str,
    params: dict,
    progress_callback: Optional[Callable] = None,
) -> List[BatchResult]:
    """Run probe design on a list of FASTA files.

    Writes output files to output_dir for each successful run.

    Args:
        fasta_paths: List of FASTA file paths to process
        output_dir: Directory for output files
        params: Dict of design parameters (passed to run_design)
        progress_callback: Optional callable(current_index, total, filename)

    Returns:
        List of BatchResult for summary table.
    """
    os.makedirs(output_dir, exist_ok=True)
    batch_results = []

    for i, fasta_path in enumerate(fasta_paths):
        fname = fasta_path.name

        if progress_callback:
            progress_callback(i, len(fasta_paths), fname)

        run_result = run_design(
            input_path=str(fasta_path),
            n_probes=params.get("n_probes", 48),
            oligo_length=params.get("oligo_length", 20),
            spacer_length=params.get("spacer_length", 2),
            target_gibbs=params.get("target_gibbs", -23.0),
            allowable_gibbs=params.get("allowable_gibbs", (-26.0, -20.0)),
            output_name=fasta_path.stem,
            species=params.get("species", "human"),
            pseudogene_mask=params.get("pseudogene_mask", False),
            genome_mask=params.get("genome_mask", False),
            index_dir=params.get("index_dir"),
            repeatmask_mode=params.get("repeatmask_mode", "none"),
            repeatmask_file=params.get("repeatmask_file"),
            save_bowtie_raw=params.get("save_bowtie_raw", False),
        )

        if run_result.result is not None and run_result.result.probes:
            prefix = fasta_path.stem
            write_output_files(
                run_result.result,
                prefix,
                mask_seqs=run_result.result.mask_strings,
                output_dir=output_dir,
            )
            batch_results.append(BatchResult(
                filename=fname,
                n_probes_found=len(run_result.result.probes),
                score=run_result.result.score,
                status="success",
                error_msg=None,
                result=run_result.result,
            ))
        else:
            batch_results.append(BatchResult(
                filename=fname,
                n_probes_found=0,
                score=float('inf'),
                status="error" if run_result.error else "no_probes",
                error_msg=run_result.error,
                result=run_result.result,
            ))

    if progress_callback:
        progress_callback(len(fasta_paths), len(fasta_paths), "Done")

    return batch_results


def format_batch_summary(results: List[BatchResult]) -> str:
    """Format batch results as a TSV summary string.

    Columns: Filename, Probes_Found, Score, Status, Error

    Args:
        results: List of BatchResult from run_batch

    Returns:
        TSV-formatted summary string
    """
    lines = ["Filename\tProbes_Found\tScore\tStatus\tError"]
    for br in results:
        score_str = f"{br.score:.4f}" if br.score != float('inf') else "N/A"
        error_str = br.error_msg or ""
        lines.append(f"{br.filename}\t{br.n_probes_found}\t{score_str}\t{br.status}\t{error_str}")
    return "\n".join(lines) + "\n"


def write_batch_summary(results: List[BatchResult], output_dir: str) -> str:
    """Write batch summary TSV to the output directory.

    Args:
        results: List of BatchResult from run_batch
        output_dir: Directory to write the summary file

    Returns:
        Path to the written summary file
    """
    summary = format_batch_summary(results)
    filepath = os.path.join(output_dir, "batch_summary.tsv")
    with open(filepath, 'w') as f:
        f.write(summary)
    return filepath
