# 🧬 ProteusIQ – Protein Intelligence & Structural Analysis

**ProteusIQ** is a comprehensive, scientifically rigorous platform for the structural and evolutionary analysis of protein sequences. It unifies biophysical property calculation, evolutionary conservation scoring, 3D structural visualization, and functional inference into a single, highly optimized web interface.

Built for robustness and scale, ProteusIQ is designed to handle complex bioinformatics pipelines seamlessly within a cloud-deployed environment.

---

## 🚀 Features

ProteusIQ executes an automated 15-step pipeline on any input sequence:

| Capability                  | Description                                                                  |
| --------------------------- | ---------------------------------------------------------------------------- |
| 📊 **Physicochemical**      | MW, pI, amino acid composition, GRAVY, instability index                     |
| 📍 **Localization**         | Signal sequence prediction, transmembrane helices, subcellular location      |
| 🔍 **Homology Search**      | Automated BLAST against the SwissProt database via NCBI + EBI APIs           |
| 🧩 **Motifs & Domains**     | PROSITE pattern scanning and InterPro domain mapping                         |
| 📈 **Conservation**         | Position-specific Shannon entropy scoring from pairwise alignments           |
| 🏗️ **3D Structure Viewer**  | Interactive 3Dmol.js viewer with 6 specialized coloring modes                |
| 🔬 **Secondary Structure**  | Experimental helix/sheet/coil determination from PDB headers                 |
| 💊 **Ligand Contacts**      | Quantitative residue-ligand interactions mapped to evolutionary conservation |
| 🎯 **Spatial Clustering**   | 3D permutation testing for conserved residue spatial clusters                |
| 🧠 **Functional Inference** | Rule-based functional prediction via an integrated reasoning engine          |
| 📐 **Multiple Alignment**   | Clustal Omega alignment generation via EBI REST APIs                         |
| 🌳 **Phylogenetics**        | Neighbor-Joining tree construction with 100-rep bootstrap support            |
| 🌀 **Intrinsic Disorder**   | TOP-IDP propensity-based unstructured region prediction                      |

---

## ⚡ Deployment Architecture & Robustness

ProteusIQ is heavily optimized for zero-configuration deployment on platforms like **Streamlit Cloud**, featuring enterprise-grade resilience to handle concurrent users and API rate limits:

- **API Rate Limiting & Resilience**: All external calls (NCBI, EBI, PDB, UniProt, AlphaFold) are routed through a centralized `resilient_request` client featuring exponential backoff, jitter, and HTTP 429 (`Retry-After`) handling.
- **Concurrency Control**: A global asyncio-compatible semaphore strictly limits outbound external API connections to 3 concurrent requests, preventing institutional IP blacklisting during heavy load.
- **Disk-Based Caching**: Integrates `diskcache` utilizing SHA-256 sequence hashing. Redundant sequence analyses bypass the 3-minute API pipeline entirely, loading instantly (0.1s) from a pre-computed local SSD cache with a 7-day TTL.
- **Zero-Footprint Memory Offloading**: Massive string payloads (like multi-megabyte PDB `.cif` files) are stripped from the Streamlit `session_state` and offloaded to persistent disk storage (`.proteusiq_cache/contents`). The frontend only reads them into memory for the exact millisecond required to render the 3D viewer, preventing RAM exhaustion crashes when users analyze dozens of proteins.
- **Automated Teardown**: Temporary PDB workspace directories are auto-cleaned immediately after the spatial clustering and ligand contact phases complete.

---

## 🏁 Quick Start

### Local Installation

ProteusIQ requires Python 3.9+.

```bash
# Clone the repository
git clone https://github.com/yourusername/proteusiq.git
cd proteusiq

# Install dependencies
pip install -r requirements.txt

# Run the application
streamlit run app.py
```

### Environment Variables

To substantially increase your NCBI BLAST rate limits (from 3 req/sec to 10 req/sec), provide an NCBI API key:

```bash
export NCBI_API_KEY="your_api_key_here"
streamlit run app.py
```

---

## 🎨 Interactive 3D Viewer

The embedded 3D viewer supports complex data overlay onto the PDB or AlphaFold structure:

| Coloring Mode           | Scientific Purpose                                                                  |
| ----------------------- | ----------------------------------------------------------------------------------- |
| **Default**             | Rainbow spectrum highlighting N-to-C terminal progression                           |
| **Conservation**        | Heatmap mapping (Blue = Variable → Red = Highly Conserved) based on Shannon entropy |
| **pLDDT Confidence**    | AlphaFold confidence mapping (Dark Blue ≥90, Light Blue ≥70, Orange ≥50, Red <50)   |
| **Ligand Contacts**     | Green highlighting for residues within 5 Å of any detected ligand                   |
| **TM Helices**          | Yellow highlighting for predicted transmembrane core segments                       |
| **Secondary Structure** | Visual distinct coloring: Red (Helix), Yellow (Sheet), Green (Coil)                 |

---

## � Scientific Methodology Summaries

- **Evolutionary Conservation**: Calculated using Shannon entropy (`H = -Σ pᵢ log₂(pᵢ)`) normalized by `log₂(20)`. The pipeline uses sample-size-adaptive Laplace smoothing on observed amino acid types. Scores < 0.2 indicate high conservation.
- **Spatial Clustering**: Identifies 3D structural clustering of conserved residues using a permutation test. It samples 1,000 random sets of equal size to calculate a p-value, determining if the top 10% most conserved residues are significantly closer in 3D space than expected by chance.
- **Phylogenetics**: Neighbor-Joining trees are constructed using a BLOSUM62 distance model, validated via 100 bootstrap replicates with majority-rule consensus.

_(Note: Secondary structure is derived strictly from experimental PDB headers; Intrinsic disorder uses the TOP-IDP scale.)_

---

## � Core Dependencies

- **Web Framework**: [Streamlit](https://streamlit.io/)
- **Bioinformatics**: [Biopython](https://biopython.org/)
- **Caching**: `diskcache`
- **Visualization**: [3Dmol.js](https://3dmol.csb.pitt.edu/), [Matplotlib](https://matplotlib.org/)
- **External Services**: NCBI BLAST, RCSB PDB, AlphaFold DB, EBI Clustal Omega, UniProt REST API


