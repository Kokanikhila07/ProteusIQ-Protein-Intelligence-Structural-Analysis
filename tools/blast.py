"""
tools/blast.py – BLAST Homology Search Module

Uses NCBI BLAST REST API (primary) and EBI BLAST REST API (fallback)
to perform blastp searches against the SwissProt database.

Strategy:
  1. Try NCBI BLAST with 5-minute polling timeout
  2. If NCBI fails/times out, fall back to EBI BLAST
  3. Parse XML results and return top 30 hits with E-value ≤ 1e-3

Scientific notes:
  - E-value is the primary significance metric, NOT identity %
  - Identity % filtering at 30% discards valid remote homologs
  - We use E-value ≤ 1e-3 as the standard cutoff (NCBI recommendation)
  - Results are sorted by E-value (ascending), not identity

Both APIs are free and require no authentication.
"""

import re
import time
import logging
import requests
from tools.api_utils import resilient_request, get_ncbi_params
import xml.etree.ElementTree as ET

logger = logging.getLogger(__name__)

# ─── NCBI Configuration ───
NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
NCBI_MAX_POLL_TIME = 300  # 5 minutes (NCBI is slow)
NCBI_POLL_INTERVAL = 10   # seconds between polls

# ─── EBI Configuration ───
EBI_BLAST_RUN = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"
EBI_BLAST_STATUS = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status"
EBI_BLAST_RESULT = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result"
EBI_MAX_POLL_TIME = 300
EBI_POLL_INTERVAL = 10

MAX_RETRIES = 2

# ─── Significance thresholds ───
MAX_EVALUE = 1e-3         # Standard BLAST significance cutoff
MIN_IDENTITY_PCT = 20.0   # Permissive identity floor (remote homologs OK)
MAX_HITS = 30             # Return up to 30 hits for MSA/conservation
NCBI_HITLIST_SIZE = 50    # Request 50 hits from NCBI to have margin after filtering


# ═════════════════════════════════════════════════════════════
# NCBI BLAST
# ═════════════════════════════════════════════════════════════

def _ncbi_submit(sequence: str) -> tuple:
    """Submit a BLAST job to NCBI. Returns (RID, RTOE)."""
    params = {
        "CMD": "Put",
        "PROGRAM": "blastp",
        "DATABASE": "swissprot",
        "QUERY": sequence,
        "FORMAT_TYPE": "XML",
        "HITLIST_SIZE": str(NCBI_HITLIST_SIZE),
        "EXPECT": str(MAX_EVALUE),
    }

    # Add NCBI identification params
    ncbi_params = get_ncbi_params()
    params.update(ncbi_params)

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            logger.info("NCBI BLAST submit (attempt %d/%d)", attempt, MAX_RETRIES)
            response = resilient_request("post", NCBI_BLAST_URL, data=params, max_retries=1)
            response.raise_for_status()
            text = response.text

            rid_match = re.search(r"RID = (\S+)", text)
            rtoe_match = re.search(r"RTOE = (\d+)", text)

            if rid_match:
                rid = rid_match.group(1)
                rtoe = int(rtoe_match.group(1)) if rtoe_match else 15
                logger.info("NCBI BLAST submitted: RID=%s, RTOE=%ds", rid, rtoe)
                return rid, rtoe

        except requests.RequestException as e:
            logger.warning("NCBI submit attempt %d failed: %s", attempt, str(e))
            if attempt < MAX_RETRIES:
                time.sleep(2 ** attempt)

    raise RuntimeError("Failed to submit NCBI BLAST job")


def _ncbi_poll(rid: str, rtoe: int, progress_callback=None) -> str:
    """Poll NCBI for results. Returns XML string."""
    initial_wait = min(rtoe, 10)
    logger.info("NCBI BLAST: waiting %ds before first poll...", initial_wait)
    time.sleep(initial_wait)

    elapsed = initial_wait

    while elapsed < NCBI_MAX_POLL_TIME:
        try:
            params = {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"}
            params.update(get_ncbi_params())
            response = resilient_request("get", NCBI_BLAST_URL, params=params, max_retries=1)
            response.raise_for_status()
            text = response.text

            if "Status=WAITING" in text:
                remaining = NCBI_MAX_POLL_TIME - elapsed
                logger.info("NCBI BLAST running (%ds elapsed, %ds remaining)...",
                           elapsed, remaining)
                if progress_callback:
                    progress_callback(f"BLAST running ({elapsed}s / {NCBI_MAX_POLL_TIME}s)...")
                time.sleep(NCBI_POLL_INTERVAL)
                elapsed += NCBI_POLL_INTERVAL
                continue
            elif "Status=FAILED" in text:
                raise RuntimeError("NCBI BLAST job failed on server")
            elif "Status=UNKNOWN" in text:
                raise RuntimeError("NCBI BLAST job expired or unknown RID")
            else:
                logger.info("NCBI BLAST results received after %ds", elapsed)
                return text

        except requests.RequestException as e:
            logger.warning("NCBI poll error: %s", str(e))
            time.sleep(NCBI_POLL_INTERVAL)
            elapsed += NCBI_POLL_INTERVAL

    raise RuntimeError(f"NCBI BLAST timed out after {NCBI_MAX_POLL_TIME}s")


def _ncbi_search(sequence: str, progress_callback=None) -> list:
    """Run NCBI BLAST end-to-end. Returns list of hits."""
    rid, rtoe = _ncbi_submit(sequence)
    xml_text = _ncbi_poll(rid, rtoe, progress_callback)
    return _parse_ncbi_xml(xml_text, len(sequence))


# ═════════════════════════════════════════════════════════════
# EBI BLAST (fallback)
# ═════════════════════════════════════════════════════════════

def _ebi_search(sequence: str, progress_callback=None) -> list:
    """
    Run BLAST via EBI REST API as fallback.
    Uses UniProtKB/Swiss-Prot database.
    """
    logger.info("Trying EBI BLAST fallback...")

    # Submit job
    try:
        params = {
            "email": "proteusiq_tool@example.com",
            "program": "blastp",
            "stype": "protein",
            "database": "uniprotkb_swissprot",
            "sequence": sequence,
            "exp": str(MAX_EVALUE),
            "alignments": str(NCBI_HITLIST_SIZE),
        }
        response = resilient_request("post", EBI_BLAST_RUN, data=params)
        response.raise_for_status()
        job_id = response.text.strip()
        logger.info("EBI BLAST submitted: %s", job_id)
    except Exception as e:
        raise RuntimeError(f"EBI BLAST submission failed: {e}")

    # Poll for results
    elapsed = 0
    time.sleep(5)  # initial wait
    elapsed = 5

    while elapsed < EBI_MAX_POLL_TIME:
        try:
            status_url = f"{EBI_BLAST_STATUS}/{job_id}"
            status_resp = resilient_request("get", status_url, timeout=(10, 15), max_retries=1)
            status = status_resp.text.strip()

            if status == "FINISHED":
                logger.info("EBI BLAST finished after %ds", elapsed)
                break
            elif status in ("FAILURE", "ERROR", "NOT_FOUND"):
                raise RuntimeError(f"EBI BLAST failed: {status}")
            else:
                logger.info("EBI BLAST %s (%ds elapsed)...", status.lower(), elapsed)
                if progress_callback:
                    progress_callback(f"EBI BLAST {status.lower()} ({elapsed}s)...")
                time.sleep(EBI_POLL_INTERVAL)
                elapsed += EBI_POLL_INTERVAL
        except requests.RequestException as e:
            logger.warning("EBI poll error: %s", e)
            time.sleep(EBI_POLL_INTERVAL)
            elapsed += EBI_POLL_INTERVAL
    else:
        raise RuntimeError(f"EBI BLAST timed out after {EBI_MAX_POLL_TIME}s")

    # Fetch XML results
    try:
        result_url = f"{EBI_BLAST_RESULT}/{job_id}/xml"
        result_resp = resilient_request("get", result_url)
        result_resp.raise_for_status()
        return _parse_ebi_xml(result_resp.text, len(sequence))
    except Exception as e:
        raise RuntimeError(f"EBI BLAST result fetch failed: {e}")


# ═════════════════════════════════════════════════════════════
# XML Parsers
# ═════════════════════════════════════════════════════════════

def _parse_ncbi_xml(xml_text: str, query_length: int) -> list:
    """Parse NCBI BLAST XML output and extract top hits."""
    hits = []

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        logger.error("Failed to parse NCBI BLAST XML: %s", str(e))
        return hits

    iterations = root.findall(".//Iteration")
    if not iterations:
        return hits

    for hit in iterations[0].findall(".//Hit"):
        hit_id = hit.findtext("Hit_id", "")
        hit_def = hit.findtext("Hit_def", "")
        hit_accession = hit.findtext("Hit_accession", "")

        organism = ""
        os_match = re.search(r"OS=(.+?)(?:\s+OX=|\s+GN=|\s+PE=|$)", hit_def)
        if os_match:
            organism = os_match.group(1).strip()

        # Parse the best HSP (High-scoring Segment Pair)
        hsp = hit.find(".//Hsp")
        if hsp is None:
            continue

        try:
            hsp_identity = int(hsp.findtext("Hsp_identity", "0"))
            hsp_positives = int(hsp.findtext("Hsp_positive", "0"))
            hsp_align_len = int(hsp.findtext("Hsp_align-len", "1"))
            hsp_evalue = float(hsp.findtext("Hsp_evalue", "999"))
            hsp_bit_score = float(hsp.findtext("Hsp_bit-score", "0"))
            hsp_query_from = int(hsp.findtext("Hsp_query-from", "1"))
            hsp_query_to = int(hsp.findtext("Hsp_query-to", "1"))

            identity_pct = round((hsp_identity / hsp_align_len) * 100, 1)
            similarity_pct = round((hsp_positives / hsp_align_len) * 100, 1)
            coverage_pct = round(
                ((hsp_query_to - hsp_query_from + 1) / query_length) * 100, 1
            )

            # Filter by E-value (primary) and identity (permissive floor)
            if hsp_evalue > MAX_EVALUE:
                continue
            if identity_pct < MIN_IDENTITY_PCT:
                continue

            hits.append({
                "hit_id": hit_id,
                "accession": hit_accession,
                "definition": hit_def,
                "organism": organism,
                "identity_pct": identity_pct,
                "similarity_pct": similarity_pct,
                "coverage_pct": coverage_pct,
                "evalue": hsp_evalue,
                "bit_score": hsp_bit_score,
            })
        except (ValueError, ZeroDivisionError):
            continue

    # Sort by E-value (best first) — scientifically correct ranking
    hits.sort(key=lambda x: x["evalue"])
    return hits[:MAX_HITS]


def _parse_ebi_xml(xml_text: str, query_length: int) -> list:
    """Parse EBI BLAST XML output and extract top hits."""
    hits = []

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        logger.error("Failed to parse EBI BLAST XML: %s", str(e))
        return hits

    # EBI uses SequenceSimilaritySearchResult format
    for hit_elem in root.iter():
        if hit_elem.tag.endswith("hit") or hit_elem.tag == "hit":
            try:
                ac = hit_elem.get("ac", "") or hit_elem.get("id", "")
                description = hit_elem.get("description", "")
                db = hit_elem.get("database", "")

                # Find alignment info
                identity_pct = 0
                evalue = 999
                align_len = 0
                bit_score = 0

                for align in hit_elem.iter():
                    if align.tag.endswith("identity") or align.tag == "identity":
                        try:
                            identity_pct = float(align.text or align.get("value", 0))
                        except (ValueError, TypeError):
                            pass
                    if align.tag.endswith("expectation") or align.tag == "expectation":
                        try:
                            evalue = float(align.text or align.get("value", 999))
                        except (ValueError, TypeError):
                            pass
                    if align.tag.endswith("alignlen") or align.tag == "alignlen":
                        try:
                            align_len = int(align.text or align.get("value", 0))
                        except (ValueError, TypeError):
                            pass
                    if align.tag.endswith("bits") or align.tag == "bits":
                        try:
                            bit_score = float(align.text or align.get("value", 0))
                        except (ValueError, TypeError):
                            pass

                # Filter by E-value and identity
                if evalue > MAX_EVALUE:
                    continue
                if identity_pct < MIN_IDENTITY_PCT or not ac:
                    continue

                coverage_pct = round(
                    (align_len / query_length) * 100, 1
                ) if align_len > 0 else 0

                # Extract organism from description
                organism = ""
                os_match = re.search(
                    r"OS=(.+?)(?:\s+OX=|\s+GN=|\s+PE=|$)", description
                )
                if os_match:
                    organism = os_match.group(1).strip()

                hits.append({
                    "hit_id": ac,
                    "accession": ac,
                    "definition": description,
                    "organism": organism,
                    "identity_pct": round(identity_pct, 1),
                    "similarity_pct": round(identity_pct, 1),  # EBI may not provide this separately
                    "coverage_pct": coverage_pct,
                    "evalue": evalue,
                    "bit_score": bit_score,
                })
            except Exception:
                continue

    # Sort by E-value (best first)
    hits.sort(key=lambda x: x["evalue"])
    return hits[:MAX_HITS]


# ═════════════════════════════════════════════════════════════
# Public API
# ═════════════════════════════════════════════════════════════

def search(sequence: str, progress_callback=None) -> dict:
    """
    Perform BLAST homology search with NCBI (primary) + EBI (fallback).

    Args:
        sequence: Uppercase single-letter amino acid string.
        progress_callback: Optional callable for progress updates.

    Returns:
        Dictionary with 'hits', 'success', 'message', 'source'.
    """
    logger.info("Starting BLAST search for sequence of length %d", len(sequence))

    # ─── Try NCBI BLAST first ───
    try:
        hits = _ncbi_search(sequence, progress_callback)
        return {
            "hits": hits,
            "success": True,
            "source": "NCBI",
            "message": (
                f"Found {len(hits)} significant hits "
                f"(E-value ≤ {MAX_EVALUE}, identity ≥ {MIN_IDENTITY_PCT}%) "
                f"via NCBI BLAST"
            ),
        }
    except RuntimeError as e:
        logger.warning("NCBI BLAST failed: %s. Trying EBI fallback...", str(e))

    # ─── Fallback to EBI BLAST ───
    try:
        hits = _ebi_search(sequence, progress_callback)
        return {
            "hits": hits,
            "success": True,
            "source": "EBI",
            "message": (
                f"Found {len(hits)} significant hits "
                f"(E-value ≤ {MAX_EVALUE}, identity ≥ {MIN_IDENTITY_PCT}%) "
                f"via EBI BLAST"
            ),
        }
    except RuntimeError as e:
        logger.error("EBI BLAST also failed: %s", str(e))
        return {
            "hits": [],
            "success": False,
            "source": None,
            "message": f"Both NCBI and EBI BLAST failed. Last error: {str(e)}",
        }
    except Exception as e:
        logger.error("Unexpected BLAST error: %s", str(e))
        return {
            "hits": [],
            "success": False,
            "source": None,
            "message": f"BLAST error: {str(e)}",
        }
