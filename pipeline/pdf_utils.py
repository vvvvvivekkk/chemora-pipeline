"""
Stage 1 of the pipeline: PDF -> page images.

Why we need this: the recognition models "look" at pictures, not at PDF files.
So the first job is to turn each page of the PDF into a plain image (PNG) that
the models can read.

We use PyMuPDF (imported as `fitz`). It's pure-Python-friendly, fast, and does
NOT need the system 'poppler' library — one less thing to break on deployment.
"""

from __future__ import annotations

import io
from dataclasses import dataclass

from PIL import Image

from . import config


@dataclass
class RenderedPage:
    """One rendered PDF page.

    index   : 0-based page number (page 0 is the first page).
    image   : the page as a PIL Image (what the models and the UI both use).
    """
    index: int
    image: Image.Image


def pdf_to_images(pdf_bytes: bytes) -> list[RenderedPage]:
    """Render a PDF (given as raw bytes) into a list of page images.

    Parameters
    ----------
    pdf_bytes : the PDF file's contents. We take bytes rather than a path so this
                works the same whether the file came from a Streamlit upload, a
                web request, or disk — the caller decides how to read it.

    Returns
    -------
    A list of RenderedPage, one per page (capped at config.MAX_PAGES).
    """
    # Imported here, not at the top, so that 'manual' mode — which may run in an
    # environment without PyMuPDF — doesn't fail at import time.
    # Newer PyMuPDF prefers `import pymupdf`; older versions only have `fitz`.
    try:
        import pymupdf as fitz
    except ImportError:
        import fitz  # older PyMuPDF (this is the same library)

    pages: list[RenderedPage] = []

    # fitz can open a PDF straight from bytes held in memory.
    with fitz.open(stream=pdf_bytes, filetype="pdf") as doc:
        total = len(doc)
        limit = total if config.MAX_PAGES == 0 else min(total, config.MAX_PAGES)

        # A "matrix" is how fitz scales the page. PDF's native resolution is 72
        # DPI, so to render at our target DPI we scale by (target / 72).
        zoom = config.PDF_RENDER_DPI / 72.0
        matrix = fitz.Matrix(zoom, zoom)

        for i in range(limit):
            page = doc[i]
            pixmap = page.get_pixmap(matrix=matrix)  # rasterize the page
            # Convert the fitz pixmap into a PIL Image via PNG bytes. PIL is what
            # the rest of the pipeline and Streamlit expect.
            png_bytes = pixmap.tobytes("png")
            image = Image.open(io.BytesIO(png_bytes)).convert("RGB")
            pages.append(RenderedPage(index=i, image=image))

    return pages
