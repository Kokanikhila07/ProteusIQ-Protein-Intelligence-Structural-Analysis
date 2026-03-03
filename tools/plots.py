"""
tools/plots.py – Visualization Plots for Streamlit

Generates charts using Streamlit's built-in charting and HTML/CSS for
custom visualizations. No extra dependencies needed.

Includes:
  - Kyte-Doolittle hydropathy profile
  - Multi-track annotated sequence viewer (scientifically rigorous)
  - Conservation heatmap bar
"""

import logging

logger = logging.getLogger(__name__)

# Kyte-Doolittle hydrophobicity scale
KD_SCALE = {
    'A':  1.8, 'C':  2.5, 'D': -3.5, 'E': -3.5, 'F':  2.8,
    'G': -0.4, 'H': -3.2, 'I':  4.5, 'K': -3.9, 'L':  3.8,
    'M':  1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
    'S': -0.8, 'T': -0.7, 'V':  4.2, 'W': -0.9, 'Y': -1.3,
}

# Amino acid property classification for coloring
AA_PROPERTIES = {
    # Hydrophobic (blue shades)
    'A': ('hydrophobic', '#5b8def'),
    'V': ('hydrophobic', '#5b8def'),
    'I': ('hydrophobic', '#5b8def'),
    'L': ('hydrophobic', '#5b8def'),
    'M': ('hydrophobic', '#5b8def'),
    'F': ('aromatic',    '#8b6bbf'),
    'W': ('aromatic',    '#8b6bbf'),
    'Y': ('aromatic',    '#8b6bbf'),
    # Polar uncharged (green shades)
    'S': ('polar', '#4CAF50'),
    'T': ('polar', '#4CAF50'),
    'N': ('polar', '#66BB6A'),
    'Q': ('polar', '#66BB6A'),
    # Positively charged (red)
    'K': ('positive', '#ef5350'),
    'R': ('positive', '#ef5350'),
    'H': ('positive', '#ef9a9a'),
    # Negatively charged (magenta)
    'D': ('negative', '#e040fb'),
    'E': ('negative', '#e040fb'),
    # Special (orange/yellow)
    'G': ('special', '#FFA726'),
    'P': ('special', '#FFB74D'),
    'C': ('special', '#FFCA28'),
}

# Three-letter codes for tooltips
THREE_LETTER = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
}


def compute_hydropathy(sequence: str, window: int = 9) -> list:
    """
    Compute Kyte-Doolittle hydropathy profile with sliding window.

    Args:
        sequence: Protein sequence.
        window: Sliding window size (default 9).

    Returns:
        List of (position, hydropathy_score) tuples.
    """
    if len(sequence) < window:
        return []

    profile = []
    half = window // 2

    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i + window]
        score = sum(KD_SCALE.get(aa, 0.0) for aa in window_seq) / window
        position = i + half + 1  # center of window, 1-based
        profile.append((position, round(score, 3)))

    return profile


def generate_sequence_annotation_html(
    sequence: str,
    entropy_scores: dict = None,
    ss_assignments: dict = None,
    tm_positions: list = None,
    conserved_positions: list = None,
    domain_annotations: list = None,
    signal_peptide: dict = None,
    disorder_scores: dict = None,
) -> str:
    """
    Generate a scientifically rigorous multi-track annotated sequence viewer.

    Tracks displayed (top to bottom):
      1. Position ruler — numbered every 10 residues
      2. Amino acid sequence — colored by physicochemical property
      3. Secondary structure — H (helix), E (sheet), C (coil) with symbols
      4. Conservation — colored gradient from conserved (red) to variable (blue)
      5. Transmembrane topology — TM helix regions highlighted
      6. Domains/motifs — colored spans for detected domains
      7. Disorder — predicted disordered regions

    Args:
        sequence: Protein sequence.
        entropy_scores: {resid: entropy} dict.
        ss_assignments: {resid: 'H'|'E'|'C'} dict.
        tm_positions: list of (start, end) tuples.
        conserved_positions: list of 1-based positions.
        domain_annotations: list of {name, start, end} dicts.
        signal_peptide: dict with 'signal_peptide' bool and optional 'cleavage_site'.
        disorder_scores: {resid: score} dict.

    Returns:
        HTML string for rendering in Streamlit.
    """
    if not sequence:
        return ""

    # Defaults
    entropy_scores = entropy_scores or {}
    ss_assignments = ss_assignments or {}
    tm_positions = tm_positions or []
    conserved_positions = conserved_positions or []
    domain_annotations = domain_annotations or []
    signal_peptide = signal_peptide or {}
    disorder_scores = disorder_scores or {}

    conserved_set = set(conserved_positions)

    # Build TM position set (convert 0-based start to 1-based)
    tm_set = set()
    for start, end in tm_positions:
        for i in range(start + 1, end + 1):
            tm_set.add(i)

    # Build domain position mapping (resid → domain_name)
    domain_map = {}
    domain_colors = [
        '#e91e63', '#9c27b0', '#673ab7', '#3f51b5', '#2196f3',
        '#009688', '#4caf50', '#ff9800', '#795548', '#607d8b',
    ]
    for idx, dom in enumerate(domain_annotations):
        color = domain_colors[idx % len(domain_colors)]
        for pos in range(dom.get("start", 0), dom.get("end", 0) + 1):
            domain_map[pos] = {
                "name": dom.get("name", dom.get("id", "Domain")),
                "color": color,
            }

    # Signal peptide range
    sp_set = set()
    if signal_peptide.get("signal_peptide"):
        cleavage = signal_peptide.get("cleavage_site", 30)
        for i in range(1, min(cleavage + 1, len(sequence) + 1)):
            sp_set.add(i)

    # Check which tracks have data
    has_ss = bool(ss_assignments)
    has_conservation = bool(entropy_scores)
    has_tm = bool(tm_set)
    has_domains = bool(domain_map)
    has_sp = bool(sp_set)
    has_disorder = bool(disorder_scores)

    lines = []
    lines.append(f"""
    <div style="font-family: 'Courier New', 'Consolas', monospace; font-size: 11px;
                line-height: 1.0; background: rgba(14,17,23,0.6);
                border: 1px solid rgba(124,157,255,0.15);
                border-radius: 10px; padding: 14px; overflow-x: auto;">
    """)

    # Process sequence in blocks of 60 residues
    block_size = 60
    for block_start in range(0, len(sequence), block_size):
        block_end = min(block_start + block_size, len(sequence))

        # ─── Track 1: Position ruler ───
        ruler_line = f'<span style="color:#4a5568;font-size:10px;">{"":>5s} </span>'
        for i in range(block_start, block_end):
            resid = i + 1
            if resid % 10 == 0:
                num_str = str(resid)
                # Right-align the number to end at this position
                ruler_line += f'<span style="color:#6b7280;font-size:10px;">{num_str[-1]}</span>'
            elif resid % 10 == 1 and resid > 1 and (resid // 10) > 0:
                # Show tens digit at position ending in 1 (for numbers that just passed)
                pass
                ruler_line += '<span style="color:#3a4550;font-size:10px;">·</span>'
            elif resid % 5 == 0:
                ruler_line += '<span style="color:#3a4550;font-size:10px;">·</span>'
            else:
                ruler_line += '<span style="font-size:10px;"> </span>'
            # Add space every 10 residues
            if (i + 1) % 10 == 0 and (i + 1) < block_end:
                ruler_line += ' '
        lines.append(ruler_line + '<br>')

        # ─── Track 2: Amino acid sequence ───
        block_num = f"{block_start + 1:5d} "
        seq_line = f'<span style="color:#6b7280;">{block_num}</span>'

        for i in range(block_start, block_end):
            resid = i + 1
            aa = sequence[i]
            prop, color = AA_PROPERTIES.get(aa, ('unknown', '#9e9e9e'))
            three = THREE_LETTER.get(aa, aa)
            hydro = KD_SCALE.get(aa, 0.0)

            # Build tooltip
            tip_parts = [f"{three}{resid} ({prop})"]
            tip_parts.append(f"Hydrophobicity: {hydro:+.1f}")
            if resid in entropy_scores:
                ent = entropy_scores[resid]
                tip_parts.append(f"Entropy: {ent:.3f}")
            if resid in conserved_set:
                tip_parts.append("★ Highly conserved")
            ss = ss_assignments.get(resid, "")
            if ss == "H":
                tip_parts.append("SS: α-helix")
            elif ss == "E":
                tip_parts.append("SS: β-sheet")
            if resid in tm_set:
                tip_parts.append("TM helix")
            if resid in domain_map:
                tip_parts.append(f"Domain: {domain_map[resid]['name']}")
            if resid in sp_set:
                tip_parts.append("Signal peptide")

            title = " | ".join(tip_parts)

            # Bold conserved residues
            weight = "bold" if resid in conserved_set else "normal"

            seq_line += (
                f'<span style="color:{color};font-weight:{weight}" '
                f'title="{title}">{aa}</span>'
            )
            if (i + 1) % 10 == 0 and (i + 1) < block_end:
                seq_line += ' '

        lines.append(seq_line + '<br>')

        # ─── Track 3: Secondary structure ───
        if has_ss:
            ss_line = f'<span style="color:#6b7280;">{"SS":>5s} </span>'
            for i in range(block_start, block_end):
                resid = i + 1
                ss = ss_assignments.get(resid, " ")
                if ss == "H":
                    char = "α"
                    color = "#ff5252"
                    bg = "rgba(255,82,82,0.15)"
                elif ss == "E":
                    char = "β"
                    color = "#ffd740"
                    bg = "rgba(255,215,64,0.15)"
                elif ss == "C":
                    char = "—"
                    color = "#66bb6a"
                    bg = "transparent"
                else:
                    char = " "
                    color = "#555"
                    bg = "transparent"
                ss_line += f'<span style="color:{color};background:{bg}">{char}</span>'
                if (i + 1) % 10 == 0 and (i + 1) < block_end:
                    ss_line += ' '
            lines.append(ss_line + '<br>')

        # ─── Track 4: Conservation ───
        if has_conservation:
            cons_line = f'<span style="color:#6b7280;">{"Cons":>5s} </span>'
            for i in range(block_start, block_end):
                resid = i + 1
                ent = entropy_scores.get(resid, None)
                if ent is not None:
                    # Gradient: conserved(red, low entropy) → variable(blue, high)
                    if ent < 0.2:
                        char = "█"
                        color = "#ff1744"
                    elif ent < 0.4:
                        char = "▓"
                        color = "#ff6e40"
                    elif ent < 0.6:
                        char = "▒"
                        color = "#ffa726"
                    elif ent < 0.8:
                        char = "░"
                        color = "#42a5f5"
                    else:
                        char = "·"
                        color = "#1e88e5"
                    cons_line += f'<span style="color:{color}" title="Entropy: {ent:.3f}">{char}</span>'
                else:
                    cons_line += '<span style="color:#333"> </span>'
                if (i + 1) % 10 == 0 and (i + 1) < block_end:
                    cons_line += ' '
            lines.append(cons_line + '<br>')

        # ─── Track 5: TM helices ───
        if has_tm or has_sp:
            topo_line = f'<span style="color:#6b7280;">{"Topo":>5s} </span>'
            for i in range(block_start, block_end):
                resid = i + 1
                if resid in sp_set:
                    topo_line += '<span style="color:#00e676;background:rgba(0,230,118,0.15)">S</span>'
                elif resid in tm_set:
                    topo_line += '<span style="color:#ffd740;background:rgba(255,215,64,0.2)">M</span>'
                else:
                    topo_line += '<span style="color:#333">·</span>'
                if (i + 1) % 10 == 0 and (i + 1) < block_end:
                    topo_line += ' '
            lines.append(topo_line + '<br>')

        # ─── Track 6: Domains ───
        if has_domains:
            dom_line = f'<span style="color:#6b7280;">{"Dom":>5s} </span>'
            for i in range(block_start, block_end):
                resid = i + 1
                if resid in domain_map:
                    d = domain_map[resid]
                    dom_line += (
                        f'<span style="color:{d["color"]};background:rgba(255,255,255,0.05)" '
                        f'title="{d["name"]}">▪</span>'
                    )
                else:
                    dom_line += '<span style="color:#333">·</span>'
                if (i + 1) % 10 == 0 and (i + 1) < block_end:
                    dom_line += ' '
            lines.append(dom_line + '<br>')

        # ─── Track 7: Disorder ───
        if has_disorder:
            dis_line = f'<span style="color:#6b7280;">{"Dis":>5s} </span>'
            for i in range(block_start, block_end):
                resid = i + 1
                dscore = disorder_scores.get(resid, None)
                if dscore is not None:
                    if dscore >= 0.5:
                        char = "~"
                        color = "#ff6f00"
                    else:
                        char = "·"
                        color = "#37474f"
                    dis_line += (
                        f'<span style="color:{color}" '
                        f'title="Disorder: {dscore:.3f}">{char}</span>'
                    )
                else:
                    dis_line += '<span style="color:#333"> </span>'
                if (i + 1) % 10 == 0 and (i + 1) < block_end:
                    dis_line += ' '
            lines.append(dis_line + '<br>')

        # Add spacing between blocks
        lines.append('<br>')

    # ─── Legend ───
    lines.append('<div style="margin-top: 6px; padding-top: 8px; border-top: 1px solid rgba(124,157,255,0.1);">')
    lines.append('<span style="color:#6b7280;font-size:10px;line-height:1.8;">')

    # Sequence colors legend
    lines.append('<b>Residues:</b> ')
    lines.append('<span style="color:#5b8def">■</span> Hydrophobic ')
    lines.append('<span style="color:#8b6bbf">■</span> Aromatic ')
    lines.append('<span style="color:#4CAF50">■</span> Polar ')
    lines.append('<span style="color:#ef5350">■</span> Positive ')
    lines.append('<span style="color:#e040fb">■</span> Negative ')
    lines.append('<span style="color:#FFA726">■</span> Special ')
    lines.append('<br>')

    if has_ss:
        lines.append('<b>SS:</b> ')
        lines.append('<span style="color:#ff5252">α</span>=helix ')
        lines.append('<span style="color:#ffd740">β</span>=sheet ')
        lines.append('<span style="color:#66bb6a">—</span>=coil ')
        lines.append('&nbsp; ')

    if has_conservation:
        lines.append('<b>Conservation:</b> ')
        lines.append('<span style="color:#ff1744">█</span>=conserved ')
        lines.append('<span style="color:#42a5f5">░</span>=moderate ')
        lines.append('<span style="color:#1e88e5">·</span>=variable ')
        lines.append('&nbsp; ')

    if has_tm or has_sp:
        lines.append('<br><b>Topology:</b> ')
        if has_sp:
            lines.append('<span style="color:#00e676">S</span>=signal peptide ')
        if has_tm:
            lines.append('<span style="color:#ffd740">M</span>=TM helix ')

    if has_disorder:
        lines.append('<b>Disorder:</b> ')
        lines.append('<span style="color:#ff6f00">~</span>=disordered ')

    lines.append('</span>')
    lines.append('</div>')
    lines.append("</div>")

    return "".join(lines)


def generate_conservation_bar_html(entropy_scores: dict, seq_length: int) -> str:
    """
    Generate a colored conservation heatmap bar.

    Each position is a thin colored segment:
      - Red = conserved (low entropy)
      - Blue = variable (high entropy)
      - Gray = no data

    Returns HTML string.
    """
    if not entropy_scores or seq_length == 0:
        return ""

    bar_width = min(seq_length, 800)
    segment_width = max(1, bar_width // seq_length)

    segments = []
    for i in range(1, seq_length + 1):
        ent = entropy_scores.get(i, 0.5)
        # Red (conserved) to blue (variable)
        r = int(255 * (1 - ent))
        b = int(255 * ent)
        g = int(60 * (1 - abs(ent - 0.5) * 2))
        color = f"rgb({r},{g},{b})"
        segments.append(
            f'<div style="display:inline-block;width:{segment_width}px;height:20px;'
            f'background:{color}" title="Pos {i}: entropy={ent:.3f}"></div>'
        )

    html = f"""
    <div style="margin: 8px 0;">
        <div style="font-size: 11px; color: #9ca3af; margin-bottom: 4px;">
            Conservation heatmap (hover for details)
        </div>
        <div style="display:flex; border-radius: 4px; overflow: hidden;
                    border: 1px solid rgba(124,157,255,0.15);">
            {"".join(segments)}
        </div>
        <div style="display: flex; justify-content: space-between;
                    font-size: 10px; color: #6b7280; margin-top: 2px;">
            <span>1</span>
            <span style="color:#ff4444">◀ Conserved</span>
            <span style="color:#4466ff">Variable ▶</span>
            <span>{seq_length}</span>
        </div>
    </div>
    """
    return html
