"""
tools/ligand.py – Ligand Detection Module

Downloads a PDB file and identifies non-water HETATM ligands bound to
the structure. Provides optional filtering of common crystallization
ions.

Uses Bio.PDB.PDBParser to parse the coordinate file.
"""

import os
import logging
import requests
import tempfile
from Bio.PDB import PDBParser

logger = logging.getLogger(__name__)

PDB_DOWNLOAD_URL = "https://files.rcsb.org/download"

# Common crystallization additives / buffer ions to optionally exclude
CRYSTALLIZATION_IONS = frozenset({
    "NA", "CL", "K", "SO4", "PO4", "GOL", "EDO", "PEG",
    "MPD", "MES", "TRS", "EPE", "ACT", "FMT", "DMS",
    "IMD", "SCN", "NO3", "BR", "IOD", "NH4",
})

# Biologically relevant metal ions (keep even when filtering)
BIOLOGICAL_METALS = frozenset({
    "ZN", "MG", "FE", "FE2", "MN", "CA", "CO", "CU",
    "CU1", "NI", "MO", "SE", "W",
})


def _download_pdb(pdb_id: str) -> str:
    """
    Download a PDB file from RCSB.

    Args:
        pdb_id: 4-character PDB identifier.

    Returns:
        Path to downloaded PDB file, or empty string on failure.
    """
    pdb_id_upper = pdb_id.upper()
    url = f"{PDB_DOWNLOAD_URL}/{pdb_id_upper}.pdb"

    try:
        logger.info("Downloading PDB file: %s", url)
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        # Save to a temporary file
        temp_dir = tempfile.mkdtemp(prefix="protein_agent_")
        filepath = os.path.join(temp_dir, f"{pdb_id_upper}.pdb")

        with open(filepath, "w") as f:
            f.write(response.text)

        logger.info("PDB file saved to %s", filepath)
        return filepath

    except requests.RequestException as e:
        logger.warning("Failed to download PDB %s: %s", pdb_id, str(e))
        return ""
    except IOError as e:
        logger.warning("Failed to save PDB file: %s", str(e))
        return ""


def detect_ligands(pdb_id: str, filter_ions: bool = False) -> dict:
    """
    Detect ligands in a PDB structure.

    Args:
        pdb_id: 4-character PDB identifier.
        filter_ions: If True, exclude common crystallization ions
                     (but keep biologically relevant metals).

    Returns:
        Dictionary with keys:
          - 'ligands': list of ligand dicts (name, chain, count)
          - 'total_unique': number of unique ligand types
          - 'filtered': whether ion filtering was applied
          - 'message': status message
    """
    logger.info("Detecting ligands in PDB %s (filter_ions=%s)", pdb_id, filter_ions)

    # Download the PDB file
    pdb_path = _download_pdb(pdb_id)
    if not pdb_path:
        return {
            "ligands": [],
            "total_unique": 0,
            "filtered": filter_ions,
            "message": f"Could not download PDB file for {pdb_id}",
        }

    try:
        # Parse the PDB file
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, pdb_path)

        # Collect HETATM residues (non-water)
        ligand_counts = {}  # (name, chain) → count

        for model in structure:
            for chain in model:
                chain_id = chain.id
                for residue in chain:
                    hetflag = residue.id[0]
                    res_name = residue.get_resname().strip()

                    # Skip standard amino acids (hetflag == ' ') and water ('W' or 'HOH')
                    if hetflag == " ":
                        continue
                    if res_name in ("HOH", "WAT", "H2O"):
                        continue

                    # Apply ion filtering if requested
                    if filter_ions:
                        if res_name in CRYSTALLIZATION_IONS and res_name not in BIOLOGICAL_METALS:
                            continue

                    key = (res_name, chain_id)
                    ligand_counts[key] = ligand_counts.get(key, 0) + 1

        # Convert to list of dicts
        ligands = []
        for (name, chain), count in sorted(ligand_counts.items()):
            ligands.append({
                "name": name,
                "chain": chain,
                "count": count,
            })

        # Get unique ligand names
        unique_names = set(lig["name"] for lig in ligands)

        logger.info("Found %d ligand entries (%d unique types) in %s",
                     len(ligands), len(unique_names), pdb_id)

        return {
            "ligands": ligands,
            "total_unique": len(unique_names),
            "filtered": filter_ions,
            "message": f"Found {len(unique_names)} unique ligand type(s) in {pdb_id}",
        }

    except Exception as e:
        logger.error("Error parsing PDB file %s: %s", pdb_id, str(e))
        return {
            "ligands": [],
            "total_unique": 0,
            "filtered": filter_ions,
            "message": f"Error parsing PDB: {str(e)}",
        }
    finally:
        # Clean up temp file
        try:
            if pdb_path and os.path.exists(pdb_path):
                os.remove(pdb_path)
                parent = os.path.dirname(pdb_path)
                if os.path.exists(parent):
                    os.rmdir(parent)
        except OSError:
            pass
