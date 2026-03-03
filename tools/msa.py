"""
tools/msa.py – Optional MSA via EBI Clustal Omega REST API

Submits sequences to the EBI Clustal Omega web service and polls for
results. Only runs when explicitly requested and at least 3 sequences
are available (query + 2 hits minimum).

Includes a clear disclaimer about external API usage.
"""

import time
import logging
import requests
from tools.api_utils import resilient_request

logger = logging.getLogger(__name__)

CLUSTALO_RUN_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
CLUSTALO_STATUS_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/status"
CLUSTALO_RESULT_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/result"

MAX_POLL_TIME = 300  # 5 minutes max
POLL_INTERVAL = 10   # seconds between polls
MIN_SEQUENCES = 3    # Clustal Omega minimum is 2, but 3 gives better results


def _build_fasta(query_sequence: str, hit_sequences: dict) -> str:
    """
    Build a multi-FASTA string from query + hit sequences.

    Args:
        query_sequence: The query protein sequence.
        hit_sequences: Dict of {accession: sequence}.

    Returns:
        Multi-FASTA formatted string.
    """
    fasta_lines = [f">Query\n{query_sequence}"]
    for acc, seq in hit_sequences.items():
        # Clean the accession for FASTA header
        clean_acc = str(acc).replace(" ", "_")[:50]
        fasta_lines.append(f">{clean_acc}\n{seq}")
    return "\n".join(fasta_lines)


def run_clustalo(
    query_sequence: str,
    hit_sequences: dict,
    email: str = "proteusiq_tool@example.com",
) -> dict:
    """
    Run Clustal Omega MSA via EBI REST API.

    Args:
        query_sequence: The query protein sequence.
        hit_sequences: Dict of {accession: sequence_string}.
        email: Email for EBI API (required by their terms).

    Returns:
        Dictionary with:
          - 'alignment': aligned FASTA string
          - 'success': bool
          - 'num_sequences': int
          - 'message': status message
          - 'disclaimer': usage note
    """
    total_seqs = 1 + len(hit_sequences)  # query + hits

    if total_seqs < MIN_SEQUENCES:
        return {
            "alignment": "",
            "success": False,
            "num_sequences": total_seqs,
            "message": f"MSA requires at least {MIN_SEQUENCES} sequences (have {total_seqs}).",
            "disclaimer": "",
        }

    fasta_input = _build_fasta(query_sequence, hit_sequences)

    logger.info("Submitting MSA job to EBI Clustal Omega (%d sequences)", total_seqs)

    # Submit job
    try:
        params = {
            "email": email,
            "sequence": fasta_input,
            "stype": "protein",
            "outfmt": "fa",
        }
        response = resilient_request("post", CLUSTALO_RUN_URL, data=params)
        response.raise_for_status()
        job_id = response.text.strip()
        logger.info("MSA job submitted: %s", job_id)
    except Exception as e:
        logger.error("MSA submission failed: %s", e)
        return {
            "alignment": "",
            "success": False,
            "num_sequences": total_seqs,
            "message": f"MSA submission failed: {e}",
            "disclaimer": "",
        }

    # Poll for results
    elapsed = 0
    while elapsed < MAX_POLL_TIME:
        time.sleep(POLL_INTERVAL)
        elapsed += POLL_INTERVAL

        try:
            status_url = f"{CLUSTALO_STATUS_URL}/{job_id}"
            status_resp = resilient_request("get", status_url, timeout=(10, 15), max_retries=1)
            status = status_resp.text.strip()

            if status == "FINISHED":
                logger.info("MSA job finished after %ds", elapsed)
                break
            elif status == "FAILURE" or status == "ERROR":
                return {
                    "alignment": "",
                    "success": False,
                    "num_sequences": total_seqs,
                    "message": f"MSA job failed with status: {status}",
                    "disclaimer": "",
                }
            elif status == "RUNNING" or status == "PENDING":
                logger.info("MSA job still %s (%ds elapsed)...", status.lower(), elapsed)
                continue
            else:
                logger.warning("Unexpected MSA status: %s", status)
                continue

        except Exception as e:
            logger.warning("MSA poll error: %s", e)
            continue
    else:
        return {
            "alignment": "",
            "success": False,
            "num_sequences": total_seqs,
            "message": f"MSA timed out after {MAX_POLL_TIME}s",
            "disclaimer": "",
        }

    # Fetch result
    try:
        result_url = f"{CLUSTALO_RESULT_URL}/{job_id}/aln-fasta"
        result_resp = resilient_request("get", result_url)
        result_resp.raise_for_status()
        alignment = result_resp.text.strip()

        logger.info("MSA alignment retrieved (%d characters)", len(alignment))

        return {
            "alignment": alignment,
            "success": True,
            "num_sequences": total_seqs,
            "message": f"MSA completed for {total_seqs} sequences",
            "disclaimer": (
                "MSA generated via Clustal Omega (EBI). "
                "May be slow for large numbers of sequences."
            ),
        }

    except Exception as e:
        logger.error("Failed to fetch MSA results: %s", e)
        return {
            "alignment": "",
            "success": False,
            "num_sequences": total_seqs,
            "message": f"Failed to retrieve MSA results: {e}",
            "disclaimer": "",
        }
