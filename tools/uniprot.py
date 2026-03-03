"""
tools/uniprot.py – UniProt Metadata Fetching

Fetches protein metadata from the UniProt REST API given a UniProt
accession. Provides gene name, protein function, subcellular location,
GO terms, and other annotations.
"""

import logging
import requests
from tools.api_utils import resilient_request

logger = logging.getLogger(__name__)

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"


def fetch_metadata(uniprot_id: str) -> dict:
    """
    Fetch protein metadata from UniProt.

    Args:
        uniprot_id: UniProt accession (e.g., 'P00533').

    Returns:
        Dictionary with gene_name, protein_name, function, subcellular_location,
        go_terms, organism, and keywords.
    """
    if not uniprot_id:
        return {"found": False, "message": "No UniProt ID provided"}

    url = f"{UNIPROT_API}/{uniprot_id}.json"

    try:
        logger.info("Fetching UniProt metadata for %s", uniprot_id)
        response = resilient_request("get", url, timeout=(10, 15))

        if response.status_code == 404:
            return {"found": False, "message": f"UniProt ID {uniprot_id} not found"}

        response.raise_for_status()
        data = response.json()

        # Gene name
        gene_name = ""
        genes = data.get("genes", [])
        if genes:
            gene_name = genes[0].get("geneName", {}).get("value", "")

        # Protein name
        protein_name = ""
        prot_desc = data.get("proteinDescription", {})
        rec_name = prot_desc.get("recommendedName", {})
        if rec_name:
            protein_name = rec_name.get("fullName", {}).get("value", "")
        elif prot_desc.get("submissionNames"):
            protein_name = prot_desc["submissionNames"][0].get("fullName", {}).get("value", "")

        # Function
        function_text = ""
        comments = data.get("comments", [])
        for comment in comments:
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    function_text = texts[0].get("value", "")
                    break

        # Subcellular location
        subcellular = []
        for comment in comments:
            if comment.get("commentType") == "SUBCELLULAR LOCATION":
                locs = comment.get("subcellularLocations", [])
                for loc in locs:
                    loc_info = loc.get("location", {})
                    loc_name = loc_info.get("value", "")
                    if loc_name:
                        subcellular.append(loc_name)

        # GO terms
        go_terms = []
        xrefs = data.get("uniProtKBCrossReferences", [])
        for xref in xrefs:
            if xref.get("database") == "GO":
                go_id = xref.get("id", "")
                props = xref.get("properties", [])
                go_term = ""
                for p in props:
                    if p.get("key") == "GoTerm":
                        go_term = p.get("value", "")
                        break
                if go_id and go_term:
                    go_terms.append({"id": go_id, "term": go_term})

        # Keywords
        keywords = [kw.get("name", "") for kw in data.get("keywords", [])]

        # Organism
        organism = data.get("organism", {}).get("scientificName", "")

        # Sequence length from UniProt
        seq_length = data.get("sequence", {}).get("length", 0)

        result = {
            "found": True,
            "gene_name": gene_name,
            "protein_name": protein_name,
            "function": function_text,
            "subcellular_location": subcellular,
            "go_terms": go_terms[:15],  # cap at 15
            "keywords": keywords[:20],
            "organism": organism,
            "sequence_length": seq_length,
        }

        logger.info("UniProt metadata fetched: %s (%s)", protein_name, gene_name)
        return result

    except requests.RequestException as e:
        logger.warning("UniProt fetch failed: %s", e)
        return {"found": False, "message": f"UniProt request failed: {e}"}
    except Exception as e:
        logger.warning("UniProt parse error: %s", e)
        return {"found": False, "message": f"Error: {e}"}
