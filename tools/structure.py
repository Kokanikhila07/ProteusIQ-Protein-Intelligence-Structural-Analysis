"""
tools/structure.py – Structure Identification Module

Searches for experimental 3D structures in the RCSB PDB using the
v2 sequence search API, with fallbacks to text search and AlphaFold DB.

Search strategy:
  1. RCSB PDB sequence search at 90% identity, e-value ≤ 1.0
  2. If no hits, retry at 30% identity, e-value ≤ 10.0
  3. If no hits, try RCSB text search by protein name from BLAST hits
  4. If no hits, try UniProt PDB cross-references
  5. If still no hits, try AlphaFold DB using UniProt ID
  6. If no UniProt ID, try UniProt ID mapping from BLAST hits

Also provides download_pdb_file() for caching PDB files for analysis.
"""

import os
import re
import logging
import tempfile
import requests
from tools.api_utils import resilient_request

logger = logging.getLogger(__name__)

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA_URL = "https://data.rcsb.org/rest/v1/core/entry"
ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction"
PDB_DOWNLOAD_URL = "https://files.rcsb.org/download"


def download_pdb_file(pdb_id: str = None, pdb_url: str = None) -> str:
    """
    Download a PDB file and return the local file path.

    Tries .cif format first (more complete), falls back to .pdb.

    Args:
        pdb_id: 4-character PDB identifier (for RCSB).
        pdb_url: Direct URL to PDB file (for AlphaFold).

    Returns:
        Path to downloaded PDB file, or empty string on failure.
    """
    if pdb_url:
        # Direct URL (e.g., AlphaFold)
        try:
            logger.info("Downloading PDB file: %s", pdb_url)
            response = resilient_request("get", pdb_url)
            response.raise_for_status()

            filename = pdb_url.split("/")[-1]
            temp_dir = tempfile.mkdtemp(prefix="proteusiq_")
            filepath = os.path.join(temp_dir, filename)

            with open(filepath, "w") as f:
                f.write(response.text)

            logger.info("PDB file saved to %s", filepath)
            return filepath
        except Exception as e:
            logger.warning("Failed to download PDB from URL: %s", e)
            return ""

    if not pdb_id:
        return ""

    pdb_id_upper = pdb_id.upper()

    # Try .pdb format first (legacy but widely supported by 3Dmol.js)
    # Try .cif format first (more complete) then .pdb
    for ext in [".cif", ".pdb"]:
        url = f"{PDB_DOWNLOAD_URL}/{pdb_id_upper}{ext}"
        logger.info("Downloading structure: %s", url)
        response = resilient_request("get", url)

        if response.status_code == 200 and len(response.text) > 100:
            temp_dir = tempfile.mkdtemp(prefix="proteusiq_")
            filename = f"{pdb_id_upper}{ext}"
            filepath = os.path.join(temp_dir, filename)

            with open(filepath, "w") as f:
                f.write(response.text)

            logger.info("Structure file saved to %s", filepath)
            return filepath
        else:
            logger.debug("Download %s%s failed with status %d", pdb_id_upper, ext, response.status_code)

    logger.warning("Failed to download structure for %s", pdb_id_upper)
    return ""


def _search_pdb_sequence(sequence: str, identity_cutoff: float = 0.9,
                         evalue_cutoff: float = 1.0) -> list:
    """
    Search RCSB PDB for structures matching the sequence.

    Args:
        sequence: Protein amino acid sequence.
        identity_cutoff: Minimum sequence identity (0.0–1.0).
        evalue_cutoff: Maximum E-value for significance.

    Returns:
        List of result dicts with pdb_id and score, or empty list.
    """
    payload = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "value": sequence,
                "identity_cutoff": identity_cutoff,
                "target": "pdb_protein_sequence",
                "evalue_cutoff": evalue_cutoff,
            },
        },
        "request_options": {
            "return_all_hits": True,
        },
        "return_type": "entry",
    }

    try:
        logger.info("Searching RCSB PDB by sequence (identity ≥%.0f%%, e-value ≤%.1f)...",
                     identity_cutoff * 100, evalue_cutoff)
        response = resilient_request("post", RCSB_SEARCH_URL, json=payload, timeout=(10, 60))

        if response.status_code == 204:
            logger.info("No PDB results at %.0f%% identity", identity_cutoff * 100)
            return []

        response.raise_for_status()
        data = response.json()

        results = []
        for result in data.get("result_set", []):
            pdb_id = result.get("identifier", "")
            score = result.get("score", 0)
            if pdb_id:
                results.append({"pdb_id": pdb_id, "score": score})

        logger.info("Found %d PDB results at %.0f%% identity",
                     len(results), identity_cutoff * 100)
        return results

    except requests.RequestException as e:
        logger.warning("PDB sequence search failed: %s", str(e))
        return []
    except (ValueError, KeyError) as e:
        logger.warning("PDB search response parse error: %s", str(e))
        return []


def _search_pdb_text(protein_name: str) -> list:
    """
    Search RCSB PDB by text (protein name) as a fallback.

    Args:
        protein_name: Protein name or description to search for.

    Returns:
        List of result dicts with pdb_id and score, or empty list.
    """
    if not protein_name or len(protein_name) < 3:
        return []

    # Clean protein name for search
    # Remove common suffixes, OS= annotations, species info
    clean_name = re.sub(r'\s*OS=.*$', '', protein_name)
    clean_name = re.sub(r'\s*\[.*?\]', '', clean_name)
    clean_name = clean_name.strip()

    if len(clean_name) < 3:
        return []

    payload = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": clean_name,
            },
        },
        "request_options": {
            "return_all_hits": False,
            "results_content_type": ["experimental"],
            "paginate": {
                "start": 0,
                "rows": 10,
            },
        },
        "return_type": "entry",
    }

    try:
        logger.info("Searching RCSB PDB by text: '%s'", clean_name[:60])
        response = resilient_request("post", RCSB_SEARCH_URL, json=payload)

        if response.status_code == 204:
            return []

        response.raise_for_status()
        data = response.json()

        results = []
        for result in data.get("result_set", []):
            pdb_id = result.get("identifier", "")
            score = result.get("score", 0)
            if pdb_id:
                results.append({"pdb_id": pdb_id, "score": score})

        logger.info("Found %d PDB results by text search", len(results))
        return results

    except Exception as e:
        logger.warning("PDB text search failed: %s", str(e))
        return []


def _search_uniprot_pdb_refs(uniprot_id: str) -> list:
    """
    Look up PDB cross-references from UniProt.

    Args:
        uniprot_id: UniProt accession.

    Returns:
        List of PDB IDs referenced by UniProt, or empty list.
    """
    if not uniprot_id:
        return []

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"

    try:
        response = resilient_request("get", url, timeout=(10, 15))
        if response.status_code != 200:
            return []

        data = response.json()
        pdb_ids = []

        for xref in data.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "PDB":
                pdb_id = xref.get("id", "")
                if pdb_id:
                    pdb_ids.append({"pdb_id": pdb_id, "score": 1.0})

        logger.info("Found %d PDB cross-references in UniProt for %s",
                     len(pdb_ids), uniprot_id)
        return pdb_ids

    except Exception as e:
        logger.debug("UniProt PDB lookup failed for %s: %s", uniprot_id, e)
        return []


def _get_pdb_details(pdb_id: str) -> dict:
    """
    Fetch experimental details for a PDB entry.

    Returns:
        Dict with method, resolution, title. Empty dict on failure.
    """
    url = f"{RCSB_DATA_URL}/{pdb_id}"

    try:
        response = resilient_request("get", url, timeout=(10, 15))
        response.raise_for_status()
        data = response.json()

        # Extract experimental method
        methods = []
        exptl = data.get("exptl", [])
        for exp in exptl:
            method = exp.get("method", "Unknown")
            methods.append(method)

        # Extract resolution
        resolution = None
        rcsb_summary = data.get("rcsb_entry_info", {})
        resolution = rcsb_summary.get("resolution_combined", [None])
        if isinstance(resolution, list) and resolution:
            resolution = resolution[0]

        # Extract title
        title = data.get("struct", {}).get("title", "")

        return {
            "method": ", ".join(methods) if methods else "Unknown",
            "resolution": resolution,
            "title": title,
        }

    except Exception as e:
        logger.warning("Failed to fetch PDB details for %s: %s", pdb_id, str(e))
        return {}


def _search_alphafold(uniprot_id: str) -> dict:
    """
    Search AlphaFold DB for a predicted structure.

    Args:
        uniprot_id: UniProt accession (e.g., 'P00533').

    Returns:
        Dict with AF model info, or empty dict on failure.
    """
    url = f"{ALPHAFOLD_API_URL}/{uniprot_id}"

    try:
        logger.info("Searching AlphaFold DB for %s...", uniprot_id)
        response = resilient_request("get", url, timeout=(10, 15))

        if response.status_code == 404:
            logger.info("No AlphaFold model for %s", uniprot_id)
            return {}

        response.raise_for_status()
        data = response.json()

        # AlphaFold API returns a list of predictions
        if isinstance(data, list) and data:
            entry = data[0]
        elif isinstance(data, dict):
            entry = data
        else:
            return {}

        pdb_url = entry.get("pdbUrl", "")
        cif_url = entry.get("cifUrl", "")

        # Get confidence score
        avg_plddt = None

        # Try the summary confidence score first
        global_metric = entry.get("globalMetricValue", None)
        if global_metric is not None:
            avg_plddt = round(float(global_metric), 1)

        # Try confidenceAvgLocalScore
        if avg_plddt is None:
            local_score = entry.get("confidenceAvgLocalScore", None)
            if local_score is not None:
                if isinstance(local_score, (float, int)):
                    avg_plddt = round(float(local_score), 1)

        # Interpret pLDDT
        interpretation = "Unknown confidence"
        if avg_plddt is not None:
            if avg_plddt >= 90:
                interpretation = "Very high confidence"
            elif avg_plddt >= 70:
                interpretation = "Confident"
            elif avg_plddt >= 50:
                interpretation = "Low confidence"
            else:
                interpretation = "Very low confidence"

        return {
            "source": "AlphaFold",
            "uniprot_id": uniprot_id,
            "pdb_url": pdb_url,
            "cif_url": cif_url,
            "avg_plddt": avg_plddt,
            "interpretation": interpretation,
        }

    except Exception as e:
        logger.warning("AlphaFold search failed for %s: %s", uniprot_id, str(e))
        return {}


def _extract_uniprot_ids(blast_hits: list) -> list:
    """
    Extract UniProt accessions from BLAST hit list.

    Returns:
        List of unique UniProt accessions.
    """
    uniprot_ids = []
    seen = set()

    for hit in (blast_hits or []):
        # Try hit_id format: sp|P12345|NAME_SPECIES
        hit_id = hit.get("hit_id", "")
        parts = hit_id.split("|")
        if len(parts) >= 2 and parts[0] in ("sp", "tr"):
            acc = parts[1]
            if acc not in seen:
                uniprot_ids.append(acc)
                seen.add(acc)
            continue

        # Try accession field
        acc = hit.get("accession", "")
        if acc and acc not in seen:
            uniprot_ids.append(acc)
            seen.add(acc)

    return uniprot_ids


def _extract_protein_name(blast_hits: list) -> str:
    """Extract a protein name from the top BLAST hit for text searching."""
    if not blast_hits:
        return ""

    definition = blast_hits[0].get("definition", "")

    # Try to extract the protein name before OS=
    name_match = re.match(r'^(.+?)\s*OS=', definition)
    if name_match:
        return name_match.group(1).strip()

    # Use first 60 chars of definition
    return definition[:60].strip()


def search(sequence: str, uniprot_id: str = None, blast_hits: list = None) -> dict:
    """
    Search for 3D structure of the protein with multi-strategy approach.

    Strategy: PDB sequence → PDB text → UniProt xrefs → AlphaFold

    Args:
        sequence: Protein amino acid sequence.
        uniprot_id: Optional UniProt accession.
        blast_hits: Optional list of BLAST hits (to extract UniProt ID and names).

    Returns:
        Dictionary with structure information, or dict with 'found': False.
    """
    logger.info("Starting structure search for sequence of length %d", len(sequence))

    # --- Step 1: Search PDB by sequence at 90% identity ---
    pdb_results = _search_pdb_sequence(sequence, identity_cutoff=0.9, evalue_cutoff=1.0)

    # --- Step 2: Retry at 30% identity with relaxed e-value ---
    if not pdb_results:
        pdb_results = _search_pdb_sequence(sequence, identity_cutoff=0.3, evalue_cutoff=10.0)

    # --- Step 3: Try UniProt PDB cross-references ---
    if not pdb_results and uniprot_id:
        pdb_results = _search_uniprot_pdb_refs(uniprot_id)

    # --- Step 3b: Try UniProt PDB cross-references from BLAST hits ---
    if not pdb_results and blast_hits:
        uniprot_ids = _extract_uniprot_ids(blast_hits)
        for uid in uniprot_ids[:3]:  # Try top 3
            pdb_results = _search_uniprot_pdb_refs(uid)
            if pdb_results:
                logger.info("Found PDB via UniProt cross-ref from BLAST hit %s", uid)
                break

    # --- Step 4: Try text search using protein name from BLAST ---
    if not pdb_results and blast_hits:
        protein_name = _extract_protein_name(blast_hits)
        if protein_name:
            pdb_results = _search_pdb_text(protein_name)

    # --- Step 5: Get details for the best PDB hit ---
    if pdb_results:
        best = pdb_results[0]
        pdb_id = best["pdb_id"]
        details = _get_pdb_details(pdb_id)

        result = {
            "found": True,
            "source": "PDB",
            "pdb_id": pdb_id,
            "method": details.get("method", "Unknown"),
            "resolution": details.get("resolution"),
            "title": details.get("title", ""),
            "score": best.get("score"),
            "total_results": len(pdb_results),
        }

        logger.info("PDB structure found: %s (%s, %.1f Å)",
                     pdb_id, result["method"],
                     result["resolution"] if result["resolution"] else 0)
        return result

    # --- Step 6: AlphaFold fallback ---
    logger.info("No PDB results, trying AlphaFold...")

    # Determine UniProt ID for AlphaFold
    af_uniprot = uniprot_id

    if not af_uniprot and blast_hits:
        uniprot_ids = _extract_uniprot_ids(blast_hits)
        if uniprot_ids:
            af_uniprot = uniprot_ids[0]
            logger.info("Using UniProt ID %s from BLAST hits for AlphaFold", af_uniprot)

    if af_uniprot:
        af_result = _search_alphafold(af_uniprot)
        if af_result:
            af_result["found"] = True
            return af_result

    # --- Step 7: Try AlphaFold with other BLAST UniProt IDs ---
    if blast_hits:
        uniprot_ids = _extract_uniprot_ids(blast_hits)
        for uid in uniprot_ids[:5]:
            if uid == af_uniprot:
                continue
            af_result = _search_alphafold(uid)
            if af_result:
                af_result["found"] = True
                logger.info("Found AlphaFold model via BLAST hit %s", uid)
                return af_result

    # --- No structure found ---
    logger.info("No structure found for this sequence")
    return {
        "found": False,
        "source": None,
        "message": "No experimental or predicted structure available.",
    }
