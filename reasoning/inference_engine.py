"""
reasoning/inference_engine.py – Functional Inference Engine

Rule-based system that aggregates evidence from all analysis modules
to infer the likely biological function of a protein. Returns a
function description, confidence level, and list of triggered rules
for transparency.

Implements 12+ rules covering major protein families:
  ABC transporters, kinases, GPCRs, secreted hormones, single-pass
  receptors, heme-binding, metalloproteins, ion channels, globins,
  proteases, DNA-binding / transcription factors, and uncharacterized.
"""

import logging

logger = logging.getLogger(__name__)


def _check_domains_contain(domains: list, keyword: str) -> bool:
    """Check if any domain name/ID contains a keyword (case-insensitive)."""
    keyword_lower = keyword.lower()
    for domain in domains:
        name = domain.get("name", "").lower()
        ps_id = domain.get("id", "").lower()
        if keyword_lower in name or keyword_lower in ps_id:
            return True
    return False


def _check_ligands_contain(ligands: list, name: str) -> bool:
    """Check if any ligand matches a given name (case-insensitive)."""
    name_upper = name.upper()
    for lig in ligands:
        if lig.get("name", "").upper() == name_upper:
            return True
    return False


def _has_metal_ion(ligands: list) -> list:
    """Return list of metal ion names found in ligands."""
    metals = {"ZN", "MG", "FE", "FE2", "MN", "CA", "CO", "CU", "CU1", "NI", "MO"}
    found = []
    for lig in ligands:
        name = lig.get("name", "").upper()
        if name in metals:
            found.append(name)
    return found


def _blast_annotation_contains(blast_hits: list, keyword: str) -> bool:
    """Check if any BLAST hit definition contains a keyword."""
    keyword_lower = keyword.lower()
    for hit in blast_hits:
        definition = hit.get("definition", "").lower()
        if keyword_lower in definition:
            return True
    return False


def _has_conserved_residue(sequence: str, conserved_positions: list, residues: str) -> bool:
    """Check if any conserved position has one of the specified residues."""
    for pos in conserved_positions:
        idx = pos - 1  # Convert 1-based to 0-based
        if 0 <= idx < len(sequence) and sequence[idx] in residues:
            return True
    return False


def infer_function(data: dict) -> dict:
    """
    Infer the biological function of a protein from aggregated analysis data.

    Args:
        data: Dictionary containing all analysis results with keys:
            - sequence (str)
            - physchem (dict)
            - tm_prediction (dict)
            - blast (dict)
            - domains (list)
            - conservation (dict)
            - structure (dict)
            - ligands (dict)

    Returns:
        Dictionary with keys:
          - function: descriptive string
          - confidence: "High", "Medium", or "Low"
          - rules_triggered: list of rule descriptions
          - evidence_count: number of supporting rules
    """
    logger.info("Running functional inference engine...")

    rules_triggered = []
    functions_inferred = []
    evidence_strength = 0
    weak_evidence_count = 0

    # Extract data components (with safe defaults)
    sequence = data.get("sequence", "")
    tm = data.get("tm_prediction", {})
    blast = data.get("blast", {})
    domains = data.get("domains", [])
    conservation = data.get("conservation", {})
    structure = data.get("structure", {})
    ligand_data = data.get("ligands", {})

    tm_helices = tm.get("tm_helices", 0)
    has_sp = tm.get("signal_peptide", False)
    blast_hits = blast.get("hits", []) if isinstance(blast, dict) else []
    ligands = ligand_data.get("ligands", []) if isinstance(ligand_data, dict) else []
    conserved_positions = conservation.get("conserved_positions", [])
    conservation_confidence = conservation.get("confidence", "Low")

    # ─────────────────────────────────────────────────────────
    # RULE 1: ABC Transporter
    # ─────────────────────────────────────────────────────────
    if tm_helices >= 6 and (
        _check_domains_contain(domains, "ABC")
        or _check_domains_contain(domains, "PS00211")
        or _blast_annotation_contains(blast_hits, "ABC")
    ):
        rules_triggered.append("Rule 1: ≥6 TM helices + ABC transporter domain/annotation")
        functions_inferred.append("Membrane ATP-binding transporter (ABC family)")
        evidence_strength += 3

    # ─────────────────────────────────────────────────────────
    # RULE 2: Kinase (ATP-binding)
    # ─────────────────────────────────────────────────────────
    if (
        _check_domains_contain(domains, "kinase")
        or _check_domains_contain(domains, "PS00107")
        or _check_domains_contain(domains, "PS00108")
    ):
        kinase_desc = "Protein kinase"
        if _check_ligands_contain(ligands, "ATP") or _check_ligands_contain(ligands, "ANP"):
            kinase_desc += " (ATP-binding confirmed)"
            evidence_strength += 3
        else:
            evidence_strength += 2

        if tm_helices == 1 and _blast_annotation_contains(blast_hits, "receptor"):
            kinase_desc = f"Receptor tyrosine kinase"
            evidence_strength += 2

        rules_triggered.append("Rule 2: Kinase domain detected")
        functions_inferred.append(kinase_desc)

    # ─────────────────────────────────────────────────────────
    # RULE 3: Secreted hormone/cytokine
    # ─────────────────────────────────────────────────────────
    if has_sp and tm_helices == 0:
        if (
            _blast_annotation_contains(blast_hits, "hormone")
            or _blast_annotation_contains(blast_hits, "cytokine")
            or _blast_annotation_contains(blast_hits, "growth factor")
        ):
            rules_triggered.append("Rule 3: Signal peptide + hormone/cytokine annotation")
            functions_inferred.append("Secreted hormone or cytokine")
            evidence_strength += 3

    # ─────────────────────────────────────────────────────────
    # RULE 4: Single-pass receptor
    # ─────────────────────────────────────────────────────────
    if tm_helices == 1 and (
        _check_domains_contain(domains, "receptor")
        or _check_domains_contain(domains, "PS00240")
        or _blast_annotation_contains(blast_hits, "receptor")
    ):
        rules_triggered.append("Rule 4: Single TM helix + receptor domain/annotation")
        functions_inferred.append("Single-pass membrane receptor")
        evidence_strength += 2

    # ─────────────────────────────────────────────────────────
    # RULE 5: Heme-binding protein
    # ─────────────────────────────────────────────────────────
    if (
        _check_ligands_contain(ligands, "HEM")
        or _check_ligands_contain(ligands, "HEC")
        or _check_domains_contain(domains, "PS00190")
    ):
        rules_triggered.append("Rule 5: Heme ligand or cytochrome domain detected")
        functions_inferred.append("Heme-binding protein")
        evidence_strength += 3

    # ─────────────────────────────────────────────────────────
    # RULE 6: Metalloprotein
    # ─────────────────────────────────────────────────────────
    metals = _has_metal_ion(ligands)
    if metals:
        metal_str = ", ".join(sorted(set(metals)))
        if _has_conserved_residue(sequence, conserved_positions, "HC"):
            rules_triggered.append(f"Rule 6: Metal ion(s) ({metal_str}) + conserved His/Cys")
            functions_inferred.append(f"Metalloprotein (binds {metal_str})")
            evidence_strength += 3
        else:
            rules_triggered.append(f"Rule 6: Metal ion(s) ({metal_str}) detected")
            functions_inferred.append(f"Metal-binding protein ({metal_str})")
            evidence_strength += 1

    # ─────────────────────────────────────────────────────────
    # RULE 7: GPCR
    # ─────────────────────────────────────────────────────────
    if tm_helices >= 7 and (
        _check_domains_contain(domains, "PS00237")
        or _blast_annotation_contains(blast_hits, "G-protein")
        or _blast_annotation_contains(blast_hits, "GPCR")
        or _blast_annotation_contains(blast_hits, "rhodopsin")
    ):
        rules_triggered.append("Rule 7: 7+ TM helices + GPCR domain/annotation")
        functions_inferred.append("G-protein coupled receptor (GPCR)")
        evidence_strength += 3

    # ─────────────────────────────────────────────────────────
    # RULE 8: Ion channel
    # ─────────────────────────────────────────────────────────
    if tm_helices >= 2 and (
        _blast_annotation_contains(blast_hits, "channel")
        or _blast_annotation_contains(blast_hits, "ion channel")
    ):
        rules_triggered.append("Rule 8: Multiple TM helices + ion channel annotation")
        functions_inferred.append("Ion channel")
        evidence_strength += 2

    # ─────────────────────────────────────────────────────────
    # RULE 9: Globin / oxygen transport
    # ─────────────────────────────────────────────────────────
    if (
        _check_domains_contain(domains, "PS01033")
        or _check_domains_contain(domains, "globin")
    ) and _check_ligands_contain(ligands, "HEM"):
        rules_triggered.append("Rule 9: Globin domain + heme ligand")
        functions_inferred.append("Globin (oxygen transport/storage)")
        evidence_strength += 3

    # ─────────────────────────────────────────────────────────
    # RULE 10: Serine protease
    # ─────────────────────────────────────────────────────────
    has_serine_motif = (
        _check_domains_contain(domains, "PS00124")
        or _check_domains_contain(domains, "PS00125")
    )
    has_protease_annotation = (
        _check_domains_contain(domains, "protease")
        or _blast_annotation_contains(blast_hits, "protease")
        or _blast_annotation_contains(blast_hits, "peptidase")
    )
    if has_serine_motif:
        rules_triggered.append("Rule 10: Serine protease motif detected")
        functions_inferred.append("Serine protease")
        evidence_strength += 3
    elif has_protease_annotation:
        rules_triggered.append("Rule 10: Protease annotation without catalytic motif")
        functions_inferred.append("Putative protease")
        evidence_strength += 1
        weak_evidence_count += 1

    # ─────────────────────────────────────────────────────────
    # RULE 11: DNA-binding / transcription factor
    # ─────────────────────────────────────────────────────────
    has_dna_binding_domain = (
        _check_domains_contain(domains, "PS00028")
        or _check_domains_contain(domains, "PS00027")
        or _check_domains_contain(domains, "zinc finger")
        or _check_domains_contain(domains, "homeobox")
    )
    if has_dna_binding_domain or _blast_annotation_contains(blast_hits, "transcription factor"):
        rules_triggered.append("Rule 11: DNA-binding domain (zinc finger/homeobox)")
        if has_dna_binding_domain:
            functions_inferred.append("DNA-binding protein / transcription factor")
            evidence_strength += 2
        else:
            functions_inferred.append("Putative DNA-binding regulator")
            evidence_strength += 1
            weak_evidence_count += 1

    # ─────────────────────────────────────────────────────────
    # RULE 12: Secreted protein (generic)
    # ─────────────────────────────────────────────────────────
    if has_sp and tm_helices == 0 and not any("Secreted" in f for f in functions_inferred):
        rules_triggered.append("Rule 12: Signal peptide, no TM helices, secreted")
        functions_inferred.append("Secreted/extracellular protein")
        evidence_strength += 1

    # ─────────────────────────────────────────────────────────
    # DEFAULT: Uncharacterized
    # ─────────────────────────────────────────────────────────
    if not functions_inferred:
        rules_triggered.append("Default: No specific rules matched")
        functions_inferred.append("Uncharacterized protein")
        evidence_strength = 0

    # ─────────────────────────────────────────────────────────
    # Determine confidence level
    # ─────────────────────────────────────────────────────────
    if evidence_strength >= 5:
        confidence = "High"
    elif evidence_strength >= 2:
        confidence = "Medium"
    else:
        confidence = "Low"

    if weak_evidence_count > 0 and confidence == "High":
        confidence = "Medium"
        rules_triggered.append("Guardrail: confidence capped due to weak annotation-only evidence")

    if conservation_confidence == "Low" and confidence == "High":
        confidence = "Medium"
        rules_triggered.append("Guardrail: downgraded due to low conservation confidence")

    # Combine all inferred functions
    function_str = "; ".join(functions_inferred)

    result = {
        "function": function_str,
        "confidence": confidence,
        "rules_triggered": rules_triggered,
        "evidence_count": len(rules_triggered),
    }

    logger.info("Functional inference: '%s' (confidence: %s, %d rules triggered)",
                function_str, confidence, len(rules_triggered))
    return result
