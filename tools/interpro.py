"""
tools/interpro.py – InterPro Domain Search Module

Queries the EBI InterPro REST API for protein domain/family annotations.
This provides much richer domain annotations than local PROSITE patterns alone,
including Pfam, SMART, CDD, PANTHER, and other integrated databases.

Strategy:
  1. Submit sequence to InterPro sequence search (InterProScan)
  2. Poll for results
  3. Parse domain hits with boundaries, families, and GO mappings

Falls back gracefully if the API is unavailable.
"""

import time
import logging
import requests
from tools.api_utils import resilient_request

logger = logging.getLogger(__name__)

IPRSCAN_RUN_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
IPRSCAN_STATUS_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status"
IPRSCAN_RESULT_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result"

MAX_POLL_TIME = 600  # InterProScan can be slow (10 min)
POLL_INTERVAL = 15   # seconds between polls


def search_domains(sequence: str, email: str = "proteusiq_tool@example.com") -> dict:
    """
    Search InterPro for domain/family annotations.

    Args:
        sequence: Protein amino acid sequence.
        email: Email for EBI API (required by terms).

    Returns:
        Dictionary with:
          - 'domains': list of domain dicts
          - 'families': list of family/superfamily dicts
          - 'go_terms': list of GO term dicts from InterPro
          - 'success': bool
          - 'message': status message
    """
    logger.info("Submitting InterProScan job for sequence of length %d", len(sequence))

    # Submit job
    try:
        params = {
            "email": email,
            "sequence": sequence,
            "stype": "p",  # protein
            "goterms": "true",
            "pathways": "true",
        }
        response = resilient_request("post", IPRSCAN_RUN_URL, data=params)
        response.raise_for_status()
        job_id = response.text.strip()
        logger.info("InterProScan job submitted: %s", job_id)
    except Exception as e:
        logger.warning("InterProScan submission failed: %s", e)
        return {
            "domains": [],
            "families": [],
            "go_terms": [],
            "success": False,
            "message": f"InterProScan submission failed: {e}",
        }

    # Poll for results
    elapsed = 0
    while elapsed < MAX_POLL_TIME:
        time.sleep(POLL_INTERVAL)
        elapsed += POLL_INTERVAL

        try:
            status_url = f"{IPRSCAN_STATUS_URL}/{job_id}"
            status_resp = resilient_request("get", status_url, timeout=(10, 15), max_retries=1)
            status = status_resp.text.strip()

            if status == "FINISHED":
                logger.info("InterProScan finished after %ds", elapsed)
                break
            elif status in ("FAILURE", "ERROR", "NOT_FOUND"):
                return {
                    "domains": [],
                    "families": [],
                    "go_terms": [],
                    "success": False,
                    "message": f"InterProScan failed: {status}",
                }
            else:
                logger.info("InterProScan %s (%ds elapsed)...", status.lower(), elapsed)
                continue
        except Exception as e:
            logger.warning("InterProScan poll error: %s", e)
            continue
    else:
        return {
            "domains": [],
            "families": [],
            "go_terms": [],
            "success": False,
            "message": f"InterProScan timed out after {MAX_POLL_TIME}s",
        }

    # Fetch JSON results
    try:
        result_url = f"{IPRSCAN_RESULT_URL}/{job_id}/json"
        result_resp = resilient_request("get", result_url)
        result_resp.raise_for_status()
        data = result_resp.json()
        return _parse_interpro_results(data)
    except Exception as e:
        logger.error("Failed to fetch InterProScan results: %s", e)
        return {
            "domains": [],
            "families": [],
            "go_terms": [],
            "success": False,
            "message": f"Failed to retrieve InterProScan results: {e}",
        }


def _parse_interpro_results(data: dict) -> dict:
    """Parse InterProScan JSON output."""
    domains = []
    families = []
    go_terms = []
    seen_go = set()

    try:
        results = data.get("results", [data]) if isinstance(data, dict) else data
        if isinstance(results, dict):
            results = [results]

        for result in results:
            matches = result.get("matches", [])

            for match in matches:
                signature = match.get("signature", {})
                sig_id = signature.get("accession", "")
                sig_name = signature.get("name", "")
                sig_desc = signature.get("description", "")

                # Source database (Pfam, SMART, Prosite, etc.)
                sig_lib = signature.get("signatureLibraryRelease", {})
                db_name = sig_lib.get("library", "Unknown")

                # InterPro entry (integrated annotation)
                entry = signature.get("entry", {})
                ipr_id = entry.get("accession", "") if entry else ""
                ipr_name = entry.get("name", "") if entry else ""
                ipr_type = entry.get("type", "") if entry else ""

                # Location(s)
                locations = match.get("locations", [])
                for loc in locations:
                    start = loc.get("start", 0)
                    end = loc.get("end", 0)
                    score = loc.get("score", None)
                    evalue = loc.get("evalue", None)

                    domain_entry = {
                        "db": db_name,
                        "signature_id": sig_id,
                        "signature_name": sig_name or sig_desc,
                        "start": start,
                        "end": end,
                        "ipr_id": ipr_id,
                        "ipr_name": ipr_name,
                    }
                    if evalue is not None:
                        domain_entry["evalue"] = evalue

                    if ipr_type in ("FAMILY", "HOMOLOGOUS_SUPERFAMILY"):
                        families.append(domain_entry)
                    else:
                        domains.append(domain_entry)

                # GO terms
                if entry and entry.get("goXRefs"):
                    for go_ref in entry["goXRefs"]:
                        go_id = go_ref.get("id", "")
                        go_name = go_ref.get("name", "")
                        go_cat = go_ref.get("category", {})
                        go_category = go_cat.get("name", "") if isinstance(go_cat, dict) else str(go_cat)

                        if go_id and go_id not in seen_go:
                            go_terms.append({
                                "id": go_id,
                                "term": go_name,
                                "category": go_category,
                            })
                            seen_go.add(go_id)

    except Exception as e:
        logger.warning("InterPro result parsing error: %s", e)

    # Sort domains by start position
    domains.sort(key=lambda x: x.get("start", 0))
    families.sort(key=lambda x: x.get("start", 0))

    total = len(domains) + len(families)
    logger.info("InterProScan found %d domains, %d families, %d GO terms",
                len(domains), len(families), len(go_terms))

    return {
        "domains": domains,
        "families": families,
        "go_terms": go_terms[:30],
        "success": True,
        "message": f"InterProScan found {total} domain/family hits",
    }
