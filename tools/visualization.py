"""
tools/visualization.py â€“ Interactive 3Dmol.js Structure Viewer

Generates an HTML snippet containing a 3Dmol.js viewer with multiple
coloring modes, toggleable representations, residue tooltips, and
interactive controls. Embedded in Streamlit via st.components.v1.html.

Coloring modes:
  - Default (element coloring)
  - Conservation (entropy gradient: blue=variable â†’ red=conserved)
  - pLDDT (AlphaFold confidence)
  - Ligand contacts (green highlight)
  - TM helices (yellow highlight)
  - Secondary structure (H=red, E=yellow, C=green)
"""

import json
import logging

logger = logging.getLogger(__name__)


def _json_safe(obj):
    """Convert numpy types to native Python types for JSON serialization."""
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    elif hasattr(obj, 'item'):  # numpy scalar (float32, int64, etc.)
        return obj.item()
    return obj


def render_advanced_viewer(
    pdb_data: str,
    pdb_format: str = "pdb",
    conservation_dict: dict = None,
    ligand_contacts: list = None,
    tm_helices: list = None,
    plddt_data: dict = None,
    secondary_structure: dict = None,
    cluster_residues: list = None,
    height: int = 600,
) -> str:
    """
    Generate an HTML snippet for an interactive 3Dmol.js protein viewer.

    Args:
        pdb_data: PDB or CIF file content as string.
        pdb_format: 'pdb' or 'cif' depending on the data format.
        conservation_dict: Optional {resid: entropy_score} for conservation coloring.
        ligand_contacts: Optional list of {resid, distance, ...} dicts.
        tm_helices: Optional list of (start, end) tuples for TM helices.
        plddt_data: Optional {resid: plddt_score} for AlphaFold coloring.
        secondary_structure: Optional {resid: 'H'|'E'|'C'} assignments.
        cluster_residues: Optional list of {resid, entropy, ...} for cluster highlighting.
        height: Viewer height in pixels.

    Returns:
        Complete HTML string for embedding via st.components.v1.html.
    """
    # Prepare JSON data for JavaScript
    js_data = {
        "conservation": conservation_dict or {},
        "ligandContacts": _prepare_ligand_contacts(ligand_contacts),
        "tmHelices": _prepare_tm_helices(tm_helices),
        "plddt": plddt_data or {},
        "secondaryStructure": secondary_structure or {},
        "clusterResidues": _prepare_cluster_residues(cluster_residues),
    }

    # Escape PDB data for JavaScript string
    pdb_escaped = json.dumps(pdb_data)
    data_json = json.dumps(_json_safe(js_data))

    html = f"""
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{ font-family: 'Inter', 'Segoe UI', sans-serif; background: transparent; }}
#viewer-container {{
    width: 100%;
    height: {height}px;
    position: relative;
    border-radius: 12px;
    overflow: hidden;
    border: 1px solid rgba(124,157,255,0.2);
    background: #0e1117;
}}
#viewer {{ width: 100%; height: 100%; }}
#controls {{
    position: absolute;
    top: 10px;
    left: 10px;
    z-index: 100;
    display: flex;
    flex-wrap: wrap;
    gap: 6px;
}}
#controls select, #controls button {{
    background: rgba(14,17,23,0.85);
    color: #e0e0e0;
    border: 1px solid rgba(124,157,255,0.3);
    border-radius: 6px;
    padding: 5px 10px;
    font-size: 12px;
    cursor: pointer;
    backdrop-filter: blur(8px);
}}
#controls select:hover, #controls button:hover {{
    border-color: rgba(124,157,255,0.6);
    background: rgba(14,17,23,0.95);
}}
#controls select:focus, #controls button:focus {{ outline: none; }}
#tooltip {{
    position: absolute;
    display: none;
    background: rgba(14,17,23,0.92);
    color: #e0e0e0;
    border: 1px solid rgba(124,157,255,0.3);
    border-radius: 8px;
    padding: 8px 12px;
    font-size: 11px;
    line-height: 1.4;
    pointer-events: none;
    z-index: 200;
    backdrop-filter: blur(8px);
    max-width: 250px;
}}
</style>
</head>
<body>
<div id="viewer-container">
    <div id="controls">
        <select id="coloring" onchange="applyColoring(this.value)" title="Coloring mode">
            <option value="default">&#x1F3A8; Default</option>
            <option value="conservation">&#x1F4C8; Conservation</option>
            <option value="plddt">&#x1F3AF; pLDDT</option>
            <option value="ligand">&#x1F49A; Ligand Contacts</option>
            <option value="tm">&#x1F7E1; TM Helices</option>
            <option value="ss">&#x1F534; Secondary Structure</option>
        </select>
        <select id="representation" onchange="applyRepresentation(this.value)" title="Representation">
            <option value="cartoon">Cartoon</option>
            <option value="surface">Surface</option>
            <option value="stick">Stick</option>
        </select>
        <button onclick="centerOnCluster()" title="Center on conserved cluster">&#x1F3AF; Cluster</button>
        <button onclick="resetView()" title="Reset view">&#x21BA; Reset</button>
    </div>
    <div id="viewer"></div>
    <div id="tooltip"></div>
</div>

<script>
(function() {{
    var pdbData = {pdb_escaped};
    var analysisData = {data_json};

    var viewer = $3Dmol.createViewer("viewer", {{
        backgroundColor: "#0e1117",
        antialias: true,
    }});

    var model = viewer.addModel(pdbData, "{pdb_format}");

    // Default style
    viewer.setStyle({{}}, {{cartoon: {{color: "spectrum"}}}});

    // Show ligands as sticks
    viewer.setStyle({{hetflag: true}}, {{stick: {{colorscheme: "greenCarbon", radius: 0.15}}}});

    viewer.zoomTo();
    viewer.render();

    // â”€â”€â”€ Tooltip on hover â”€â”€â”€
    var tooltip = document.getElementById("tooltip");
    viewer.setHoverable({{}}, true,
        function(atom, viewer, event, container) {{
            if (!atom) return;
            var resi = atom.resi;
            var resn = atom.resn;
            var chain = atom.chain;
            var info = "<b>" + resn + " " + resi + "</b> (Chain " + chain + ")";

            // Conservation info
            var cons = analysisData.conservation;
            if (cons && cons[resi] !== undefined) {{
                var entropy = cons[resi];
                var cls = entropy < 0.2 ? "Highly conserved" :
                          entropy < 0.5 ? "Moderate" : "Variable";
                info += "<br>Entropy: " + entropy.toFixed(3) + " (" + cls + ")";
            }}

            // Ligand contact info
            var lc = analysisData.ligandContacts;
            var contactKey = chain + ":" + resi;
            if (lc && lc[contactKey] !== undefined) {{
                info += "<br>Ligand contact: " + lc[contactKey].toFixed(1) + " Ã…";
            }} else if (lc && lc[String(resi)] !== undefined) {{
                info += "<br>Ligand contact: " + lc[String(resi)].toFixed(1) + " Ã…";
            }}

            tooltip.innerHTML = info;
            tooltip.style.display = "block";
            tooltip.style.left = (event.offsetX + 15) + "px";
            tooltip.style.top = (event.offsetY + 15) + "px";
        }},
        function(atom) {{
            tooltip.style.display = "none";
        }}
    );

    // â”€â”€â”€ Coloring functions â”€â”€â”€
    window.applyColoring = function(mode) {{
        // Reset protein style first
        viewer.setStyle({{hetflag: false}}, {{}});

        switch(mode) {{
            case "conservation":
                applyConservationColoring();
                break;
            case "plddt":
                applyPlddtColoring();
                break;
            case "ligand":
                applyLigandColoring();
                break;
            case "tm":
                applyTmColoring();
                break;
            case "ss":
                applySsColoring();
                break;
            default:
                viewer.setStyle({{hetflag: false}}, {{cartoon: {{color: "spectrum"}}}});
        }}

        // Always show ligands
        viewer.setStyle({{hetflag: true}}, {{stick: {{colorscheme: "greenCarbon", radius: 0.15}}}});
        viewer.render();
    }};

    function applyConservationColoring() {{
        var cons = analysisData.conservation;
        if (!cons || Object.keys(cons).length === 0) {{
            viewer.setStyle({{hetflag: false}}, {{cartoon: {{color: "spectrum"}}}});
            return;
        }}

        // Default to gray
        viewer.setStyle({{hetflag: false}}, {{cartoon: {{color: "#888888"}}}});

        for (var resi in cons) {{
            var entropy = cons[resi];
            // Gradient: conserved (red, low entropy) â†’ variable (blue, high entropy)
            var r = Math.round(255 * (1 - entropy));
            var b = Math.round(255 * entropy);
            var g = Math.round(80 * (1 - Math.abs(entropy - 0.5) * 2));
            var color = "rgb(" + r + "," + g + "," + b + ")";
            viewer.setStyle({{resi: parseInt(resi), hetflag: false}},
                            {{cartoon: {{color: color}}}});
        }}
    }}

    function applyPlddtColoring() {{
        var plddt = analysisData.plddt;
        if (!plddt || Object.keys(plddt).length === 0) {{
            viewer.setStyle({{hetflag: false}}, {{cartoon: {{color: "spectrum"}}}});
            return;
        }}

        viewer.setStyle({{hetflag: false}}, {{cartoon: {{color: "#888888"}}}});

        for (var resi in plddt) {{
            var score = plddt[resi];
            var color;
            if (score >= 90) color = "#0053D6";       // dark blue
            else if (score >= 70) color = "#65CBF3";   // light blue
            else if (score >= 50) color = "#FFDB13";   // orange
            else color = "#FF7D45";                     // red

            viewer.setStyle({{resi: parseInt(resi), hetflag: false}},
                            {{cartoon: {{color: color}}}});
        }}
    }}

    function applyLigandColoring() {{
        var lc = analysisData.ligandContacts;
        // Default style
        viewer.setStyle({{hetflag: false}}, {{cartoon: {{color: "#555555"}}}});

        if (lc && Object.keys(lc).length > 0) {{
            for (var key in lc) {{
                var parts = key.split(":");
                var selector = null;
                if (parts.length === 2) {{
                    selector = {{chain: parts[0], resi: parseInt(parts[1]), hetflag: false}};
                }} else {{
                    selector = {{resi: parseInt(key), hetflag: false}};
                }}
                viewer.setStyle(selector,
                                {{cartoon: {{color: "#00ff88"}}, stick: {{color: "#00ff88", radius: 0.1}}}});
            }}
        }}
    }}

    function applyTmColoring() {{
        viewer.setStyle({{hetflag: false}}, {{cartoon: {{color: "#555555"}}}});
        var tmh = analysisData.tmHelices;
        if (tmh && tmh.length > 0) {{
            for (var i = 0; i < tmh.length; i++) {{
                var start = tmh[i][0];
                var end = tmh[i][1];
                for (var r = start; r <= end; r++) {{
                    viewer.setStyle({{resi: r, hetflag: false}},
                                    {{cartoon: {{color: "#FFD700"}}}});
                }}
            }}
        }}
    }}

    function applySsColoring() {{
        var ss = analysisData.secondaryStructure;
        if (!ss || Object.keys(ss).length === 0) {{
            viewer.setStyle({{hetflag: false}}, {{cartoon: {{color: "spectrum"}}}});
            return;
        }}

        viewer.setStyle({{hetflag: false}}, {{cartoon: {{color: "#555555"}}}});

        for (var resi in ss) {{
            var type = ss[resi];
            var color;
            if (type === "H") color = "#FF4444";       // helix = red
            else if (type === "E") color = "#FFD700";   // sheet = yellow
            else color = "#44CC44";                      // coil = green

            viewer.setStyle({{resi: parseInt(resi), hetflag: false}},
                            {{cartoon: {{color: color}}}});
        }}
    }}

    // â”€â”€â”€ Representation â”€â”€â”€
    window.applyRepresentation = function(rep) {{
        var currentColoring = document.getElementById("coloring").value;

        switch(rep) {{
            case "surface":
                viewer.setStyle({{hetflag: false}}, {{cartoon: {{opacity: 0.5, color: "spectrum"}}}});
                viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity: 0.7, color: "white"}}, {{hetflag: false}});
                break;
            case "stick":
                viewer.removeAllSurfaces();
                viewer.setStyle({{hetflag: false}}, {{stick: {{colorscheme: "Jmol", radius: 0.1}}}});
                break;
            default:
                viewer.removeAllSurfaces();
                applyColoring(currentColoring);
                return;
        }}

        viewer.setStyle({{hetflag: true}}, {{stick: {{colorscheme: "greenCarbon", radius: 0.15}}}});
        viewer.render();
    }};

    // â”€â”€â”€ Center on cluster â”€â”€â”€
    window.centerOnCluster = function() {{
        var clusterResidues = analysisData.clusterResidues;
        if (clusterResidues && clusterResidues.length > 0) {{
            var first = clusterResidues[0];
            var zoomSel = first.chain
                ? {{chain: first.chain, resi: first.resid}}
                : {{resi: first.resid}};
            viewer.zoomTo(zoomSel);

            // Highlight cluster as spheres
            for (var i = 0; i < clusterResidues.length; i++) {{
                var c = clusterResidues[i];
                var sel = c.chain
                    ? {{chain: c.chain, resi: c.resid, hetflag: false}}
                    : {{resi: c.resid, hetflag: false}};
                viewer.addStyle(sel,
                               {{sphere: {{color: "#FF4488", radius: 0.6, opacity: 0.6}}}});
            }}
            viewer.render();
        }}
    }};

    // â”€â”€â”€ Reset â”€â”€â”€
    window.resetView = function() {{
        viewer.removeAllSurfaces();
        document.getElementById("coloring").value = "default";
        document.getElementById("representation").value = "cartoon";
        viewer.setStyle({{}}, {{cartoon: {{color: "spectrum"}}}});
        viewer.setStyle({{hetflag: true}}, {{stick: {{colorscheme: "greenCarbon", radius: 0.15}}}});
        viewer.zoomTo();
        viewer.render();
    }};

}})();
</script>
</body>
</html>
"""

    return html


def _prepare_ligand_contacts(contacts: list) -> dict:
    """Convert ligand contact list to {chain:resid: min_distance} for JS."""
    if not contacts:
        return {}

    result = {}
    for contact in contacts:
        resid = contact.get("resid")
        chain = contact.get("chain", "")
        dist = contact.get("distance", 999)
        if resid is not None:
            key = f"{chain}:{resid}" if chain else str(resid)
            if key not in result or dist < result[key]:
                result[key] = dist
    return result


def _prepare_cluster_residues(cluster_residues: list) -> list:
    """Convert cluster residue list to [{chain, resid}, ...] for JS."""
    if not cluster_residues:
        return []

    prepared = []
    for residue in cluster_residues:
        resid = residue.get("resid")
        if resid is None:
            continue
        prepared.append({
            "chain": residue.get("chain", ""),
            "resid": resid,
        })
    return prepared


def _prepare_tm_helices(helices: list) -> list:
    """Convert TM helix positions to [[start, end], ...] for JS."""
    if not helices:
        return []
    return [[s + 1, e] for s, e in helices]  # convert to 1-based

