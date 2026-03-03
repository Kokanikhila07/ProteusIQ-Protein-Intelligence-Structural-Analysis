"""
tools/disorder.py – Intrinsic Disorder Prediction

Predicts intrinsically disordered regions (IDRs) using amino acid
composition-based propensity scoring. This is a simplified IUPred-like
approach based on the TOP-IDP scale (Campen et al., 2008).

Method:
  1. Assign per-residue disorder propensity from the TOP-IDP scale
  2. Apply sliding window smoothing (default: 21 residues)
  3. Classify regions with score > 0.5 as disordered

The TOP-IDP scale ranks amino acids by their propensity to promote
intrinsic disorder, derived from statistical analysis of disordered
protein databases (DisProt).

Reference:
  Campen et al. (2008) "TOP-IDP-Scale: A New Amino Acid Scale
  Measuring Propensity for Intrinsic Disorder" Protein Pept Lett 15:956-963
"""

import logging

logger = logging.getLogger(__name__)

# TOP-IDP disorder propensity scale
# Higher values → more disorder-promoting
# Normalized to approximately [0, 1] range
TOP_IDP_SCALE = {
    'A':  0.06, 'R':  0.18, 'N':  0.23, 'D':  0.19, 'C': -0.02,
    'E':  0.74, 'Q':  0.56, 'G':  0.17, 'H': -0.01, 'I': -0.49,
    'L': -0.34, 'K':  0.59, 'M': -0.10, 'F': -0.42, 'P':  0.99,
    'S':  0.34, 'T':  0.05, 'W': -0.49, 'Y': -0.24, 'V': -0.29,
}

DISORDER_THRESHOLD = 0.5  # Score above this → disordered
DEFAULT_WINDOW = 21       # Standard IUPred window size


def predict(sequence: str, window: int = DEFAULT_WINDOW) -> dict:
    """
    Predict intrinsically disordered regions.

    Args:
        sequence: Protein amino acid sequence.
        window: Sliding window size (default 21).

    Returns:
        Dictionary with:
          - 'scores': dict {resid_1based: score}
          - 'disordered_regions': list of (start, end) tuples (1-based)
          - 'disorder_content': float, fraction of disordered residues
          - 'num_disordered': int
          - 'summary': human-readable summary
    """
    seq_len = len(sequence)
    if seq_len < 10:
        return {
            "scores": {},
            "disordered_regions": [],
            "disorder_content": 0.0,
            "num_disordered": 0,
            "summary": "Sequence too short for disorder prediction.",
        }

    # Raw per-residue propensity
    raw_scores = []
    for aa in sequence:
        raw_scores.append(TOP_IDP_SCALE.get(aa, 0.0))

    # Sliding window smoothing
    half_w = window // 2
    smoothed = {}

    for i in range(seq_len):
        start = max(0, i - half_w)
        end = min(seq_len, i + half_w + 1)
        window_scores = raw_scores[start:end]
        avg = sum(window_scores) / len(window_scores)
        # Normalize to [0, 1] range using sigmoid-like transformation
        normalized = _normalize_score(avg)
        smoothed[i + 1] = round(normalized, 4)  # 1-based

    # Find disordered regions (contiguous runs above threshold)
    disordered_regions = []
    in_disorder = False
    region_start = 0

    for i in range(1, seq_len + 1):
        if smoothed[i] >= DISORDER_THRESHOLD:
            if not in_disorder:
                region_start = i
                in_disorder = True
        else:
            if in_disorder:
                if i - region_start >= 5:  # minimum 5 residues for a region
                    disordered_regions.append((region_start, i - 1))
                in_disorder = False

    # Handle case where disorder extends to the end
    if in_disorder and (seq_len + 1 - region_start) >= 5:
        disordered_regions.append((region_start, seq_len))

    # Statistics
    n_disordered = sum(1 for s in smoothed.values() if s >= DISORDER_THRESHOLD)
    disorder_content = round(n_disordered / seq_len * 100, 1)

    # Summary
    if disordered_regions:
        region_strs = [f"{s}-{e}" for s, e in disordered_regions]
        summary = (
            f"{disorder_content}% predicted disordered ({n_disordered}/{seq_len} residues). "
            f"{len(disordered_regions)} disordered region(s): {', '.join(region_strs)}"
        )
    else:
        summary = f"Protein appears largely ordered ({disorder_content}% disordered content)."

    logger.info("Disorder prediction: %s", summary)

    return {
        "scores": smoothed,
        "disordered_regions": disordered_regions,
        "disorder_content": disorder_content,
        "num_disordered": n_disordered,
        "summary": summary,
    }


def _normalize_score(raw_score: float) -> float:
    """
    Normalize raw disorder score to [0, 1] range.
    Uses a linear mapping from the approximate TOP-IDP range [-0.5, 1.0].
    """
    # Map [-0.5, 1.0] → [0.0, 1.0]
    normalized = (raw_score + 0.5) / 1.5
    return max(0.0, min(1.0, normalized))
