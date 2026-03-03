"""
tools/alignment.py â€“ MSA & Conservation Module

Computes per-position conservation scores from pairwise alignments between
the query sequence and top BLAST hits. Uses Biopython's PairwiseAligner
with BLOSUM62 substitution matrix.

Conservation is measured as the fraction of aligned hits that have an
**identical** residue at each query position. Similarity-based scoring
is handled downstream by the Shannon entropy module.

Sequence Fetching Strategy:
  1. Try UniProt REST API (primary â€” BLAST hits come from SwissProt)
  2. Fall back to NCBI efetch if UniProt fails
"""

import re
import logging
import requests
import time
from tools.api_utils import resilient_request
from Bio.Align import PairwiseAligner, substitution_matrices

logger = logging.getLogger(__name__)

MAX_FETCH_HITS = 20   # Fetch sequences for top 20 hits
MIN_SEQS_FOR_CONSERVATION = 5  # Minimum unique homologs for meaningful conservation
HIGH_CONFIDENCE_MIN_SEQS = 10


def _extract_uniprot_accession(hit: dict) -> str:
    """
    Extract a UniProt accession from a BLAST hit.

    Handles multiple formats:
      - sp|P00533|EGFR_HUMAN  â†’ P00533
      - tr|A0A0G2JMB2|...     â†’ A0A0G2JMB2
      - Direct accession field â†’ as-is
    """
    # Try accession field first
    accession = hit.get("accession", "")

    # Try extracting from hit_id (sp|XXXXX|NAME format)
    if not accession or len(accession) < 4:
        hit_id = hit.get("hit_id", "")
        parts = hit_id.split("|")
        if len(parts) >= 2 and parts[0] in ("sp", "tr"):
            accession = parts[1]

    return accession.strip()


def _fetch_sequence_uniprot(accession: str) -> str:
    """
    Fetch a protein sequence from UniProt REST API.

    This is the preferred method since BLAST hits come from SwissProt.

    Args:
        accession: UniProt accession (e.g., 'P00533').

    Returns:
        Protein sequence string, or empty string on failure.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"

    try:
        response = resilient_request("get", url, headers={"Accept": "text/plain"})
        if response.status_code == 200:
            lines = response.text.strip().split("\n")
            seq_lines = [line.strip() for line in lines if not line.startswith(">")]
            sequence = "".join(seq_lines).upper()
            if len(sequence) > 10:
                logger.info("Fetched sequence from UniProt for %s (length %d)",
                           accession, len(sequence))
                return sequence
        logger.debug("UniProt fetch for %s returned status %d", accession, response.status_code)
    except Exception as e:
        logger.debug("UniProt fetch failed for %s: %s", accession, str(e))

    return ""


def _fetch_sequence_ncbi(accession: str) -> str:
    """
    Fetch a protein sequence from NCBI by accession (fallback).

    Args:
        accession: NCBI protein accession (e.g., 'P00533').

    Returns:
        Protein sequence string, or empty string on failure.
    """
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "protein",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    }

    try:
        response = resilient_request("get", url, params=params, timeout=(10, 20))
        response.raise_for_status()
        lines = response.text.strip().split("\n")
        seq_lines = [line.strip() for line in lines if not line.startswith(">")]
        sequence = "".join(seq_lines).upper()
        if len(sequence) > 10:
            logger.info("Fetched sequence from NCBI for %s (length %d)",
                       accession, len(sequence))
            return sequence
    except Exception as e:
        logger.debug("NCBI fetch failed for %s: %s", accession, str(e))

    return ""


def _fetch_sequence(accession: str) -> str:
    """
    Fetch a protein sequence, trying UniProt first, then NCBI.

    Args:
        accession: Protein accession.

    Returns:
        Protein sequence string, or empty string on failure.
    """
    # Try UniProt first (BLAST hits from SwissProt â†’ UniProt is authoritative)
    seq = _fetch_sequence_uniprot(accession)
    if seq:
        return seq

    # Fall back to NCBI
    seq = _fetch_sequence_ncbi(accession)
    return seq


def _align_sequences(query: str, subject: str) -> dict:
    """
    Perform global pairwise alignment and return position mapping.

    Returns:
        Dict mapping query position (0-based) â†’ subject residue character.
    """
    try:
        aligner = PairwiseAligner()
        aligner.mode = "global"

        # Load BLOSUM62 substitution matrix
        try:
            aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        except Exception:
            logger.warning("Could not load BLOSUM62, using default scoring")
            aligner.match_score = 2
            aligner.mismatch_score = -1

        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5

        # Get the best alignment
        alignments = aligner.align(query, subject)
        if not alignments:
            return {}

        best = alignments[0]

        # Use the alignment coordinates to map positions
        aligned = best.aligned

        # aligned is a tuple of two arrays, each with shape (n_blocks, 2)
        query_blocks = aligned[0]
        subject_blocks = aligned[1]

        # Build a mapping of query positions to aligned subject positions
        query_to_subject = {}
        for (q_start, q_end), (s_start, s_end) in zip(query_blocks, subject_blocks):
            for offset in range(q_end - q_start):
                q_pos = q_start + offset
                s_pos = s_start + offset
                if q_pos < len(query) and s_pos < len(subject):
                    query_to_subject[q_pos] = subject[s_pos]

        return query_to_subject

    except Exception as e:
        logger.warning("Pairwise alignment failed: %s", str(e))
        return {}


def compute_conservation(query_sequence: str, blast_hits: list) -> dict:
    """
    Compute per-position conservation from BLAST hit pairwise alignments.

    Args:
        query_sequence: The query protein sequence.
        blast_hits: List of BLAST hit dictionaries (from blast.py).

    Returns:
        Dictionary with keys:
          - 'conserved_positions': list of 1-based positions with conservation >= 0.8
          - 'conservation_scores': list of (position, score) for all positions
          - 'position_alignments': dict {pos: [aa_list]} for entropy calculation
          - 'hit_sequences': dict {accession: sequence} for MSA
          - 'summary': human-readable summary string
          - 'num_sequences': number of sequences used
          - 'skipped': bool, True if conservation was skipped
          - 'message': status message
    """
    logger.info(
        "Computing conservation for query of length %d with %d BLAST hits",
        len(query_sequence),
        len(blast_hits),
    )

    if len(blast_hits) < MIN_SEQS_FOR_CONSERVATION:
        return {
            "conserved_positions": [],
            "conservation_scores": [],
            "position_alignments": {},
            "hit_sequences": {},
            "summary": (
                f"Conservation analysis skipped: fewer than "
                f"{MIN_SEQS_FOR_CONSERVATION} homologous sequences available."
            ),
            "num_sequences": len(blast_hits),
            "confidence": "Low",
            "quality_flags": ["insufficient_homologs"],
            "skipped": True,
            "message": "Insufficient hits for conservation analysis",
        }

    top_hits = blast_hits[:MAX_FETCH_HITS]
    hit_sequences = {}
    seen_sequences = {query_sequence}
    fetch_failures = 0

    for hit in top_hits:
        accession = _extract_uniprot_accession(hit)
        if not accession:
            fetch_failures += 1
            continue

        if accession in hit_sequences:
            continue

        seq = _fetch_sequence(accession)
        if seq:
            if seq not in seen_sequences:
                hit_sequences[accession] = seq
                seen_sequences.add(seq)
        else:
            fetch_failures += 1

        time.sleep(0.3)

    logger.info(
        "Fetched %d/%d unique hit sequences (%d failures)",
        len(hit_sequences),
        len(top_hits),
        fetch_failures,
    )

    if len(hit_sequences) < MIN_SEQS_FOR_CONSERVATION:
        return {
            "conserved_positions": [],
            "conservation_scores": [],
            "position_alignments": {},
            "hit_sequences": hit_sequences,
            "summary": (
                f"Conservation analysis skipped: could only fetch "
                f"{len(hit_sequences)} unique homolog sequences "
                f"(need >= {MIN_SEQS_FOR_CONSERVATION})."
            ),
            "num_sequences": len(hit_sequences),
            "confidence": "Low",
            "quality_flags": ["insufficient_unique_homologs"],
            "skipped": True,
            "message": "Could not fetch sufficient sequences",
        }

    position_matches = [0] * len(query_sequence)
    all_mappings = []

    for i, (acc, hit_seq) in enumerate(hit_sequences.items()):
        logger.info("Aligning hit %d/%d (%s)...", i + 1, len(hit_sequences), acc)
        mapping = _align_sequences(query_sequence, hit_seq)
        all_mappings.append(mapping)

        if isinstance(mapping, dict):
            for q_pos, s_aa in mapping.items():
                if s_aa == query_sequence[q_pos]:
                    position_matches[q_pos] += 1

    position_alignments = {}
    for pos in range(len(query_sequence)):
        aas = [query_sequence[pos]]
        for mapping in all_mappings:
            if isinstance(mapping, dict) and pos in mapping:
                aas.append(mapping[pos])
        position_alignments[pos] = aas

    n_seqs = len(hit_sequences)
    conservation_scores = []
    conserved_positions = []

    for pos in range(len(query_sequence)):
        score = round(position_matches[pos] / n_seqs, 2)
        conservation_scores.append((pos + 1, score))
        if score >= 0.8:
            conserved_positions.append(pos + 1)

    total = len(query_sequence)
    n_conserved = len(conserved_positions)
    pct = round(n_conserved / total * 100, 1) if total > 0 else 0

    if n_seqs >= HIGH_CONFIDENCE_MIN_SEQS:
        confidence = "High"
    elif n_seqs >= MIN_SEQS_FOR_CONSERVATION:
        confidence = "Medium"
    else:
        confidence = "Low"

    quality_flags = []
    if n_seqs < HIGH_CONFIDENCE_MIN_SEQS:
        quality_flags.append("limited_homolog_depth")
    if fetch_failures > 0:
        quality_flags.append("partial_sequence_fetch_failures")

    if n_conserved > 20:
        pos_str = (
            ", ".join(str(p) for p in conserved_positions[:15])
            + f", ... ({n_conserved} total)"
        )
    else:
        pos_str = ", ".join(str(p) for p in conserved_positions) if conserved_positions else "none"

    summary = (
        f"{n_conserved} positions ({pct}%) are highly conserved (>=80% identity) "
        f"across {n_seqs} unique homologs. Confidence: {confidence}. "
        f"Conserved positions: {pos_str}"
    )

    logger.info(
        "Conservation analysis complete: %d/%d positions conserved across %d homologs",
        n_conserved,
        total,
        n_seqs,
    )

    return {
        "conserved_positions": conserved_positions,
        "conservation_scores": conservation_scores,
        "position_alignments": position_alignments,
        "hit_sequences": hit_sequences,
        "summary": summary,
        "num_sequences": n_seqs,
        "confidence": confidence,
        "quality_flags": quality_flags,
        "skipped": False,
        "message": "Conservation analysis completed successfully",
    }

