"""
tools/motif.py – Motif & Domain Detection Module

Scans a protein sequence against a curated set of PROSITE patterns.
Uses a PROSITE-to-regex converter that handles all standard PROSITE
pattern syntax elements.

PROSITE pattern syntax conversion:
  x    → .       (any amino acid)
  (n)  → {n}     (repeat n times)
  (n,m)→ {n,m}   (repeat n to m times)
  [ABC]→ [ABC]   (one of)
  {ABC}→ [^ABC]  (none of)
  <    → ^       (N-terminus)
  >    → $       (C-terminus)
  -    → (separator, removed)

This module ships with a curated set of ~25 important PROSITE patterns.
It can be extended to parse the full prosite.dat file from ExPASy.
"""

import re
import os
import logging

logger = logging.getLogger(__name__)


def prosite_to_regex(pattern: str) -> str:
    """
    Convert a PROSITE pattern string to a Python regular expression.

    Args:
        pattern: PROSITE pattern string (e.g., '[RK]-x(2,3)-[DE]-x-C').

    Returns:
        Python regex string.

    Examples:
        >>> prosite_to_regex('[RK]-x(2,3)-[DE]-x-C')
        '[RK].{2,3}[DE].C'
    """
    # Remove the trailing period if present (PROSITE convention)
    pattern = pattern.rstrip(".")

    # Replace the PROSITE separator '-'
    # First, handle elements that use hyphens internally (shouldn't exist, but be safe)
    # We process token by token

    result = ""
    i = 0
    tokens = pattern.split("-")

    for token in tokens:
        if not token:
            continue

        # Handle N-terminus anchor
        if token == "<":
            result += "^"
            continue
        # Handle C-terminus anchor
        if token == ">":
            result += "$"
            continue

        # Handle anchors attached to tokens
        if token.startswith("<"):
            result += "^"
            token = token[1:]
        ends_with_anchor = False
        if token.endswith(">"):
            ends_with_anchor = True
            token = token[:-1]

        # Handle 'x' (any amino acid)
        if token.lower() == "x":
            result += "."
        # Handle 'x(n)' or 'x(n,m)' – any AA repeated
        elif token.lower().startswith("x(") and token.endswith(")"):
            repeat = token[2:-1]
            result += ".{" + repeat + "}"
        # Handle [ABC] – character class (keep as-is)
        elif token.startswith("[") and "]" in token:
            bracket_end = token.index("]") + 1
            char_class = token[:bracket_end]
            remainder = token[bracket_end:]
            result += char_class
            # Handle optional repeat like [ABC](2)
            if remainder.startswith("(") and remainder.endswith(")"):
                repeat = remainder[1:-1]
                result += "{" + repeat + "}"
        # Handle {ABC} – negated character class
        elif token.startswith("{") and "}" in token:
            brace_end = token.index("}") + 1
            chars = token[1:brace_end - 1]
            remainder = token[brace_end:]
            result += "[^" + chars + "]"
            if remainder.startswith("(") and remainder.endswith(")"):
                repeat = remainder[1:-1]
                result += "{" + repeat + "}"
        # Handle single amino acid letter
        elif len(token) == 1 and token.isalpha():
            result += token.upper()
        # Handle single letter with repeat like A(3)
        elif len(token) >= 3 and token[0].isalpha() and token[1] == "(" and token.endswith(")"):
            letter = token[0].upper()
            repeat = token[2:-1]
            result += letter + "{" + repeat + "}"
        else:
            # Fallback: use as-is (shouldn't happen for valid PROSITE patterns)
            result += token

        if ends_with_anchor:
            result += "$"

    return result


# ─────────────────────────────────────────────────────────────
# Curated PROSITE patterns for common protein domains/motifs
# Each entry: (PROSITE_ID, name, PROSITE_pattern)
# ─────────────────────────────────────────────────────────────
CURATED_PATTERNS = [
    # Kinase & phosphorylation
    ("PS00107", "Protein kinase C-terminal domain",
     "[LIVMFYC]-x-[HY]-x-D-[LIVMFY]-K-x(2)-N-[LIVMFYCT](3)"),
    ("PS00108", "Protein kinase ATP-binding region",
     "[LIVMFYC]-G-{P}-G-{P}-[FYWMGSTNH]-[SGA]-{PW}-[LIVCAT]-{PD}-x-[GSTACLIVMFY]-x(5,18)-[LIVMFYWCSTAR]-[APTS]-[LIVMFYSTAGCQ]-K"),
    ("PS00001", "N-glycosylation site", "N-{P}-[ST]-{P}"),
    ("PS00004", "cAMP/cGMP-dependent protein kinase phosphorylation site",
     "[RK](2)-x-[ST]"),
    ("PS00005", "Protein kinase C phosphorylation site", "[ST]-x-[RK]"),
    ("PS00006", "Casein kinase II phosphorylation site", "[ST]-x(2)-[DE]"),
    ("PS00008", "N-myristoylation site", "G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}"),
    # Receptor motifs
    ("PS00240", "Receptor L domain signature",
     "C-x(5)-C-x(3)-[LIVMFYAH]-x(3)-[LIVMFYA]-x(3,4)-C"),
    ("PS50311", "EGF-like domain", "C-x(2,7)-C-x(3,14)-C-x(1,3)-[GSTAPIMVQH]-x(0,6)-C"),
    # Zinc finger
    ("PS00028", "Zinc finger C2H2 type",
     "C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H"),
    # Enzyme active sites
    ("PS00124", "Serine protease active site (His)",
     "[LIVMFYAG]-[LIVMSTAG]-[LIVFYMS]-D-[ST]-G-[STAG]"),
    ("PS00125", "Serine protease active site (Ser)",
     "G-D-S-G-G-[LIVMFY]"),
    # Metal binding
    ("PS00196", "EF-hand calcium-binding domain",
     "D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]"),
    # Structural motifs
    ("PS00010", "Leucine zipper pattern",
     "L-x(6)-L-x(6)-L-x(6)-L"),
    ("PS00015", "Nuclear localization signal",
     "K-R-[LIVMFYAT]-x-K-[RK]"),
    # Transport
    ("PS00211", "ABC transporter signature",
     "[LIVMFY]-S-[GS]-G-x(3)-[RKA]-[LIVMYA]-x-[LIVMF]-[AG]"),
    # Immunoglobulin
    ("PS00290", "Immunoglobulin domain",
     "[FY]-x-C-x-[VA]-x-H"),
    # GPCR
    ("PS00237", "G-protein coupled receptor signature",
     "[GSTALIVMYWC]-[GSTANCPDE]-{EDPKRH}-x(2)-[LIVMNQGA]-x(2)-[LIVMFT]-[GSTANC]-[LIVMFYWSTAC]-[DENH]-R-[FYWCSH]-x(2)-[LIVM]"),
    # Cytochrome / heme
    ("PS00190", "Cytochrome c family heme-binding site",
     "C-{CPWHF}-{CPWR}-C-H-{CFYW}"),
    # Thioredoxin
    ("PS00194", "Thioredoxin family active site",
     "[LIVMF]-[LIVMSTAGY]-x-[LIVMFYC]-x-x-C-x-x-C-[LIVMFYA]"),
    # AAA ATPase
    ("PS00674", "ATP synthase alpha/beta family",
     "[LIVMFYAG]-G-x(4)-[LIVMFA]-G-D-[LIVMGA]-x-[GASC]-[LIVMFAT]-x-R"),
    # Homeobox
    ("PS00027", "Homeobox domain",
     "[LIVMFYG]-{SH}-[LIVMFYGS]-[EQ]-[LIVMA]-x-x-[LIVMFYA]-x(3)-[LIVMF]-x(2)-K-[LIVMF]-x(3)-[LIVMF]"),
    # Globin
    ("PS01033", "Globin family profile",
     "[LIVMFYW]-x-[LIVMFYW]-x(3)-H-x(3)-[LIVMFY]-x(3,5)-[LIVMFYW]"),
    # WD repeat
    ("PS00678", "WD repeat signature",
     "[LIVMSTAC]-x(2)-[LIVMFYWSTAGC]-[FYWSTAGC]-[DENSHQ]-x-[LIVMSTAG]-W-D"),
    # SH2 domain
    ("PS50001", "SH2 domain",
     "[LIVMFY]-x-[LIVMFY]-x(3)-[LIVMFYC]-x(2)-[GSTALIV]-x(5,10)-[LIVMFYW]-x-[LIVMFYW]-x-[LIVMFYW]"),
]


def _compile_patterns() -> list:
    """
    Compile curated PROSITE patterns into regex objects.

    Returns:
        List of (id, name, compiled_regex) tuples.
    """
    compiled = []
    for ps_id, name, pattern in CURATED_PATTERNS:
        try:
            regex_str = prosite_to_regex(pattern)
            regex_obj = re.compile(regex_str)
            compiled.append((ps_id, name, regex_obj))
        except re.error as e:
            logger.warning("Failed to compile pattern %s (%s): %s", ps_id, name, str(e))
    return compiled


# Pre-compile all patterns at module load time
_COMPILED_PATTERNS = _compile_patterns()


def scan(sequence: str) -> list:
    """
    Scan a protein sequence for PROSITE motif/domain matches.

    Args:
        sequence: Uppercase single-letter amino acid string.

    Returns:
        List of match dictionaries with keys: id, name, start, end.
        Positions are 1-based.
    """
    logger.info("Scanning sequence (length %d) against %d PROSITE patterns",
                len(sequence), len(_COMPILED_PATTERNS))

    matches = []

    for ps_id, name, regex in _COMPILED_PATTERNS:
        for match in regex.finditer(sequence):
            matches.append({
                "id": ps_id,
                "name": name,
                "start": match.start() + 1,  # Convert to 1-based
                "end": match.end(),           # 1-based inclusive
            })

    # Sort by start position
    matches.sort(key=lambda x: x["start"])

    logger.info("Found %d motif/domain matches", len(matches))
    return matches
