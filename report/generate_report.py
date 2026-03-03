"""
report/generate_report.py – HTML Report Generator

Renders the Jinja2 template with analysis data to produce a
downloadable HTML report. Optional PDF generation via WeasyPrint
(falls back to HTML if WeasyPrint is unavailable).
"""

import os
import logging
from jinja2 import Environment, FileSystemLoader

logger = logging.getLogger(__name__)

# Path to the template directory (same directory as this file)
TEMPLATE_DIR = os.path.dirname(os.path.abspath(__file__))


def generate_html(data: dict) -> str:
    """
    Render the protein analysis report as an HTML string.

    Args:
        data: Dictionary containing all analysis results. Expected keys:
            - sequence_name (str)
            - physchem (dict)
            - tm_prediction (dict)
            - blast (dict)
            - domains (list)
            - conservation (dict)
            - structure (dict)
            - ligands (dict)  → mapped to 'ligands_data' in template
            - inference (dict)

    Returns:
        Rendered HTML string.
    """
    logger.info("Generating HTML report...")

    try:
        env = Environment(
            loader=FileSystemLoader(TEMPLATE_DIR),
            autoescape=True,
        )
        template = env.get_template("template.html")

        # Map data keys to template variable names
        context = {
            "sequence_name": data.get("sequence_name", "Query Protein"),
            "physchem": data.get("physchem", {}),
            "tm_prediction": data.get("tm_prediction", {}),
            "blast": data.get("blast", {}),
            "domains": data.get("domains", []),
            "conservation": data.get("conservation", {}),
            "structure": data.get("structure", {}),
            "ligands_data": data.get("ligands", {}),
            "inference": data.get("inference", {}),
        }

        html = template.render(**context)
        logger.info("HTML report generated successfully (%d chars)", len(html))
        return html

    except Exception as e:
        logger.error("Failed to generate HTML report: %s", str(e))
        raise


def generate_pdf(data: dict) -> bytes:
    """
    Generate a PDF report from the analysis data.

    Uses WeasyPrint if available; raises ImportError otherwise.

    Args:
        data: Same dictionary as generate_html.

    Returns:
        PDF content as bytes.

    Raises:
        ImportError: If WeasyPrint is not installed.
    """
    try:
        from weasyprint import HTML as WeasyHTML
    except ImportError:
        logger.warning("WeasyPrint not available – PDF generation disabled")
        raise ImportError(
            "WeasyPrint is not installed. Install with: pip install weasyprint "
            "(requires system dependencies). Falling back to HTML download."
        )

    html_string = generate_html(data)
    logger.info("Generating PDF report with WeasyPrint...")
    pdf_bytes = WeasyHTML(string=html_string).write_pdf()
    logger.info("PDF report generated (%d bytes)", len(pdf_bytes))
    return pdf_bytes
