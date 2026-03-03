"""
tools/physchem.py – Physicochemical Analysis Module

Computes physicochemical properties of a protein sequence using
Biopython's ProtParam (Bio.SeqUtils.ProtParam).

Returns: length, molecular weight (Da), theoretical pI, amino acid
composition (counts and percentages), aromaticity, instability index,
and GRAVY (grand average of hydropathy).
"""

import logging
from Bio.SeqUtils.ProtParam import ProteinAnalysis

logger = logging.getLogger(__name__)


def analyze(sequence: str) -> dict:
    """
    Compute physicochemical properties of a protein sequence.

    Args:
        sequence: Uppercase single-letter amino acid string.

    Returns:
        Dictionary with keys: length, molecular_weight, theoretical_pI,
        amino_acid_counts, amino_acid_percent, aromaticity,
        instability_index, gravy.
    """
    logger.info("Starting physicochemical analysis for sequence of length %d", len(sequence))

    try:
        analysis = ProteinAnalysis(sequence)

        # Amino acid counts and percentages
        aa_counts = analysis.count_amino_acids()
        aa_percent = analysis.get_amino_acids_percent()
        # Round percentages to 4 decimal places for readability
        aa_percent_rounded = {aa: round(pct, 4) for aa, pct in aa_percent.items()}

        result = {
            "length": len(sequence),
            "molecular_weight": round(analysis.molecular_weight(), 2),
            "theoretical_pI": round(analysis.isoelectric_point(), 2),
            "amino_acid_counts": aa_counts,
            "amino_acid_percent": aa_percent_rounded,
            "aromaticity": round(analysis.aromaticity(), 4),
            "instability_index": round(analysis.instability_index(), 2),
            "gravy": round(analysis.gravy(), 4),
        }

        logger.info(
            "Physicochemical analysis complete: MW=%.2f Da, pI=%.2f",
            result["molecular_weight"],
            result["theoretical_pI"],
        )
        return result

    except Exception as e:
        logger.error("Physicochemical analysis failed: %s", str(e))
        raise
