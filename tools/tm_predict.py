"""
tools/tm_predict.py – Signal Peptide & Transmembrane Prediction Module

Uses Kyte-Doolittle hydrophobicity scale with sliding window analysis to
predict signal peptides and transmembrane helices. Derives subcellular
localization from the combination of these predictions.

Signal peptide heuristic:
  - Scan first 70 residues with window=10, KD threshold=1.5
  - Look for cleavage site motif (small neutral residues A,G,S,C,T at -1,-3)
  
Transmembrane helix prediction:
  - Sliding window of 19, KD threshold=1.5, merge overlapping helices
"""

import re
import logging

logger = logging.getLogger(__name__)

# Kyte-Doolittle hydrophobicity scale
KD_SCALE = {
    'A':  1.8, 'C':  2.5, 'D': -3.5, 'E': -3.5, 'F':  2.8,
    'G': -0.4, 'H': -3.2, 'I':  4.5, 'K': -3.9, 'L':  3.8,
    'M':  1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
    'S': -0.8, 'T': -0.7, 'V':  4.2, 'W': -0.9, 'Y': -1.3,
}


def _window_hydrophobicity(sequence: str, start: int, window_size: int) -> float:
    """Compute average KD hydrophobicity for a window of residues."""
    window = sequence[start:start + window_size]
    total = sum(KD_SCALE.get(aa, 0.0) for aa in window)
    return total / len(window)


def _predict_signal_peptide(sequence: str) -> dict:
    """
    Predict signal peptide in the first 70 residues.

    Returns:
        Dict with 'has_signal_peptide' (bool), 'cleavage_site' (int or None),
        and 'hydrophobic_region' (tuple or None).
    """
    result = {
        "has_signal_peptide": False,
        "cleavage_site": None,
        "hydrophobic_region": None,
    }

    # Only scan the first 70 residues (or entire seq if shorter)
    scan_region = sequence[:min(70, len(sequence))]
    window_size = 10

    if len(scan_region) < window_size:
        return result

    # Step 1: Find hydrophobic region using sliding window
    hydrophobic_start = None
    hydrophobic_end = None

    for i in range(len(scan_region) - window_size + 1):
        avg_kd = _window_hydrophobicity(scan_region, i, window_size)
        if avg_kd > 1.5:
            if hydrophobic_start is None:
                hydrophobic_start = i
            hydrophobic_end = i + window_size
        elif hydrophobic_start is not None:
            # End of hydrophobic stretch
            break

    if hydrophobic_start is None:
        return result

    result["hydrophobic_region"] = (hydrophobic_start, hydrophobic_end)

    # Step 2: Look for cleavage site motif after hydrophobic region
    # Small neutral residues (A, G, S, C, T) at positions -1 and -3 
    # relative to the cleavage site
    search_start = hydrophobic_end
    search_end = min(search_start + 30, len(sequence))
    search_region = sequence[search_start:search_end]

    # Pattern: position -3 is [AGSCT], position -2 is any, position -1 is [AGSCT]
    cleavage_pattern = re.compile(r'[AGSCT].[AGSCT]')
    match = cleavage_pattern.search(search_region)

    if match:
        # Cleavage site is after the matched triplet
        cleavage_pos = search_start + match.end()
        result["has_signal_peptide"] = True
        result["cleavage_site"] = cleavage_pos
        logger.info(
            "Signal peptide predicted: hydrophobic region %d-%d, cleavage at %d",
            hydrophobic_start, hydrophobic_end, cleavage_pos,
        )
    else:
        logger.info("Hydrophobic region found but no cleavage site motif detected")

    return result


def _predict_tm_helices(sequence: str) -> list:
    """
    Predict transmembrane helices using sliding window of 19 and KD scale.

    Returns:
        List of (start, end) tuples for predicted TM helices (merged).
    """
    window_size = 19
    threshold = 1.5

    if len(sequence) < window_size:
        return []

    # Find all windows exceeding the threshold
    tm_windows = []
    for i in range(len(sequence) - window_size + 1):
        avg_kd = _window_hydrophobicity(sequence, i, window_size)
        if avg_kd > threshold:
            tm_windows.append((i, i + window_size))

    if not tm_windows:
        return []

    # Merge overlapping windows
    merged = [tm_windows[0]]
    for start, end in tm_windows[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end:
            # Overlapping – extend the previous helix
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))

    logger.info("Predicted %d transmembrane helix/helices", len(merged))
    return merged


def predict(sequence: str) -> dict:
    """
    Predict signal peptide, transmembrane helices, and subcellular localization.

    Args:
        sequence: Uppercase single-letter amino acid string.

    Returns:
        Dictionary with keys: signal_peptide (bool), cleavage_site (int|None),
        tm_helices (int), tm_helix_positions (list of tuples),
        localization (str).
    """
    logger.info("Starting TM/signal peptide prediction for sequence of length %d", len(sequence))

    sp_result = _predict_signal_peptide(sequence)
    tm_positions = _predict_tm_helices(sequence)

    has_sp = sp_result["has_signal_peptide"]
    num_tm = len(tm_positions)

    # Determine localization based on signal peptide + TM helix combination
    if has_sp and num_tm == 0:
        localization = "Secreted"
    elif has_sp and num_tm > 0:
        localization = "Membrane (secretory pathway)"
    elif not has_sp and num_tm > 0:
        localization = "Membrane (internal)"
    else:
        localization = "Cytosolic"

    result = {
        "signal_peptide": has_sp,
        "cleavage_site": sp_result["cleavage_site"],
        "hydrophobic_region": sp_result["hydrophobic_region"],
        "tm_helices": num_tm,
        "tm_helix_positions": tm_positions,
        "localization": localization,
    }

    logger.info(
        "TM prediction complete: signal_peptide=%s, tm_helices=%d, localization=%s",
        has_sp, num_tm, localization,
    )
    return result
