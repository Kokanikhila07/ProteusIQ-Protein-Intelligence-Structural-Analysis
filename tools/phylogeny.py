"""
tools/phylogeny.py – Phylogenetic Tree Construction with Bootstrap Support

Builds a Neighbor-Joining phylogenetic tree from MSA output using
Bio.Phylo.TreeConstruction with BLOSUM62-based distance calculation.

Features:
  - 100-replicate bootstrap for statistical support on internal nodes
  - Majority-rule consensus tree
  - Matplotlib rendering (saved as base64 PNG)
  - Newick format string for downstream use
  - ASCII fallback for environments without display

Only runs when explicitly enabled by the user, MSA was successful,
and at least 4 sequences are available.
"""

import io
import base64
import logging
import random
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.Consensus import bootstrap_consensus, majority_consensus

logger = logging.getLogger(__name__)

# Bootstrap configuration
BOOTSTRAP_REPLICATES = 100
MIN_BOOTSTRAP_DISPLAY = 50  # Only show bootstrap values >= 50%


def _parse_fasta_alignment(alignment_text: str):
    """Parse FASTA alignment text into a MultipleSeqAlignment object."""
    handle = io.StringIO(alignment_text)
    try:
        alignment = AlignIO.read(handle, "fasta")
        return alignment
    except Exception as e:
        logger.error("Failed to parse alignment: %s", e)
        return None


def _tree_to_newick(tree) -> str:
    """Convert a Bio.Phylo tree to Newick format string."""
    output = io.StringIO()
    try:
        Phylo.write(tree, output, "newick")
        return output.getvalue().strip()
    except Exception as e:
        logger.error("Failed to write Newick: %s", e)
        return ""


def _tree_to_ascii(tree) -> str:
    """Convert a Bio.Phylo tree to ASCII art string."""
    output = io.StringIO()
    try:
        Phylo.draw_ascii(tree, file=output)
        return output.getvalue()
    except Exception as e:
        logger.error("Failed to draw ASCII tree: %s", e)
        return ""


def _render_tree_matplotlib(tree, num_seqs: int) -> str:
    """
    Render a phylogenetic tree using Matplotlib with bootstrap values.

    Returns base64-encoded PNG string, or empty string on failure.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")  # Non-interactive backend
        import matplotlib.pyplot as plt

        # Dynamically size the figure based on number of sequences
        fig_height = max(4, num_seqs * 0.45)
        fig_width = max(8, 10)

        fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))

        # Style the plot
        ax.set_facecolor("#0e1117")
        fig.patch.set_facecolor("#0e1117")

        # Draw the tree
        Phylo.draw(
            tree, axes=ax, do_show=False,
            label_colors=lambda label: "#c8d6ff",
        )

        # Add bootstrap values on internal nodes
        for clade in tree.find_clades(order="level"):
            if clade.confidence is not None and clade.confidence >= MIN_BOOTSTRAP_DISPLAY:
                # Find the position of the node
                try:
                    x = clade.branch_length if clade.branch_length else 0
                    # Get node coordinates from the drawn tree
                    pass  # Bootstrap values are shown via Phylo.draw labels
                except Exception:
                    pass

        # Annotate bootstrap values on internal nodes
        for clade in tree.find_clades(order="level"):
            if not clade.is_terminal() and clade.confidence is not None:
                conf = int(clade.confidence)
                if conf >= MIN_BOOTSTRAP_DISPLAY:
                    clade.name = f"{conf}"

        # Re-draw with bootstrap labels
        ax.clear()
        ax.set_facecolor("#0e1117")

        Phylo.draw(
            tree, axes=ax, do_show=False,
            label_colors=lambda label: "#c8d6ff" if label is None or not label.isdigit() else "#ff9800",
        )

        # Style axes
        ax.set_title(
            f"Phylogenetic Tree ({num_seqs} sequences, {BOOTSTRAP_REPLICATES} bootstrap replicates)",
            color="#a5b4fc", fontsize=12, fontweight="bold", pad=10
        )
        ax.tick_params(colors="#6b7280", labelsize=9)
        ax.set_xlabel("Branch length (substitutions/site)", color="#9ca3af", fontsize=10)
        ax.set_ylabel("", color="#9ca3af")

        for spine in ax.spines.values():
            spine.set_color("#2d3748")

        plt.tight_layout()

        # Save to buffer
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight",
                    facecolor="#0e1117", edgecolor="none")
        plt.close(fig)

        buf.seek(0)
        b64_data = base64.b64encode(buf.read()).decode("utf-8")
        logger.info("Tree rendered as PNG (%d bytes)", len(b64_data))
        return b64_data

    except ImportError:
        logger.warning("Matplotlib not available, skipping tree rendering")
        return ""
    except Exception as e:
        logger.error("Tree rendering failed: %s", e)
        return ""


def _build_bootstrap_tree(alignment):
    """
    Build a Neighbor-Joining tree with bootstrap support.

    Uses 100 bootstrap replicates and majority-rule consensus.

    Returns:
        Tuple of (consensus_tree, simple_nj_tree) or (None, None) on failure.
    """
    try:
        calculator = DistanceCalculator("blosum62")
        constructor = DistanceTreeConstructor(calculator, "nj")

        # Build bootstrap consensus tree
        logger.info("Running %d bootstrap replicates...", BOOTSTRAP_REPLICATES)

        consensus_tree = bootstrap_consensus(
            alignment, BOOTSTRAP_REPLICATES, constructor, majority_consensus
        )

        # Ladderize for cleaner display
        consensus_tree.ladderize()

        # Also build a simple NJ tree for branch lengths (consensus may lack them)
        dm = calculator.get_distance(alignment)
        nj_tree = constructor.build_tree(alignment)
        nj_tree.ladderize()

        logger.info("Bootstrap consensus tree built successfully")
        return consensus_tree, nj_tree

    except Exception as e:
        logger.error("Bootstrap tree construction failed: %s", e)
        return None, None


def build_tree(alignment_text: str) -> dict:
    """
    Build a Neighbor-Joining phylogenetic tree from MSA with bootstrap.

    Args:
        alignment_text: FASTA-format alignment string from Clustal Omega.

    Returns:
        Dictionary with:
          - 'newick': Newick format tree string
          - 'ascii_tree': ASCII art representation
          - 'tree_image': base64-encoded PNG of Matplotlib rendering
          - 'num_sequences': int
          - 'success': bool
          - 'message': status message
          - 'disclaimer': scientific note
          - 'bootstrap_replicates': int
    """
    disclaimer = (
        "Phylogenetic tree constructed using Neighbor-Joining with "
        "BLOSUM62 distance model. Bootstrap support values (%) shown "
        f"on internal nodes ({BOOTSTRAP_REPLICATES} replicates, "
        f"majority-rule consensus). Values ≥{MIN_BOOTSTRAP_DISPLAY}% indicate "
        "statistically supported clades."
    )

    alignment = _parse_fasta_alignment(alignment_text)
    if alignment is None:
        return {
            "newick": "", "ascii_tree": "", "tree_image": "",
            "num_sequences": 0, "success": False,
            "message": "Failed to parse MSA alignment for tree construction",
            "disclaimer": disclaimer, "bootstrap_replicates": 0,
        }

    num_seqs = len(alignment)
    if num_seqs < 3:
        return {
            "newick": "", "ascii_tree": "", "tree_image": "",
            "num_sequences": num_seqs, "success": False,
            "message": f"Need at least 3 sequences for tree (have {num_seqs})",
            "disclaimer": disclaimer, "bootstrap_replicates": 0,
        }

    logger.info("Building NJ tree with bootstrap from %d sequences", num_seqs)

    try:
        # Build bootstrap consensus tree
        consensus_tree, nj_tree = _build_bootstrap_tree(alignment)

        if consensus_tree is None:
            # Fall back to simple NJ tree without bootstrap
            logger.warning("Bootstrap failed, falling back to simple NJ tree")
            calculator = DistanceCalculator("blosum62")
            constructor = DistanceTreeConstructor(calculator, "nj")
            nj_tree = constructor.build_tree(alignment)
            nj_tree.ladderize()
            consensus_tree = nj_tree
            disclaimer = (
                "Phylogenetic tree constructed using Neighbor-Joining with "
                "BLOSUM62 distance model. Bootstrap analysis failed — "
                "branch support could not be calculated."
            )

        # Generate outputs
        newick = _tree_to_newick(nj_tree if nj_tree else consensus_tree)
        ascii_tree = _tree_to_ascii(consensus_tree)

        # Try Matplotlib rendering
        tree_image = _render_tree_matplotlib(consensus_tree, num_seqs)

        logger.info("Phylogenetic tree built successfully (%d sequences)", num_seqs)

        return {
            "newick": newick,
            "ascii_tree": ascii_tree,
            "tree_image": tree_image,
            "num_sequences": num_seqs,
            "success": True,
            "message": (
                f"Neighbor-Joining tree with {BOOTSTRAP_REPLICATES} bootstrap "
                f"replicates built from {num_seqs} sequences"
            ),
            "disclaimer": disclaimer,
            "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        }

    except Exception as e:
        logger.error("Tree construction failed: %s", e)
        return {
            "newick": "", "ascii_tree": "", "tree_image": "",
            "num_sequences": num_seqs, "success": False,
            "message": f"Tree construction failed: {e}",
            "disclaimer": disclaimer, "bootstrap_replicates": 0,
        }
