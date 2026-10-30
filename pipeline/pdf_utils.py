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
from dataclasses import dataclass, field

from PIL import Image

from . import config


@dataclass
class RenderedPage:
    """One rendered PDF page.

    index   : 0-based page number (page 0 is the first page).
    image   : the page as a PIL Image (what the models and the UI both use).
    plus_signs   : normalised (x, y) centres of '+' signs drawn as vector lines
                   (charges in schemes, which the models ignore).
    label_images : normalised (x1, y1, x2, y2) boxes of glyph-sized pictures
                   embedded in the page. Some publishers paste R-group labels
                   (R', R'' in a coloured circle/box) into vector drawings as
                   tiny images; post-processing uses these boxes to find them.
    """
    index: int
    image: Image.Image
    label_images: list = field(default_factory=list)
    plus_signs: list = field(default_factory=list)  # (x, y) of '+' drawn as two crossed lines


# An embedded image smaller than this (fraction of page width and height) is a
# label glyph, not a figure: a whole structure drawn as an image is far bigger.
LABEL_IMAGE_MAX = 0.05


def _plus_signs(page) -> list[tuple[float, float]]:
    """Centres of '+' signs drawn as two stroked lines: one horizontal, one
    vertical, equal length, crossing at their midpoints, glyph-sized."""
    import math
    W, H = page.rect.width, page.rect.height
    found = []
    for d in page.get_drawings():
        items = d["items"]
        if d.get("fill") is not None or len(items) != 2 or any(it[0] != "l" for it in items):
            continue
        (a1, a2), (b1, b2) = [(it[1], it[2]) for it in items]
        la, lb = math.dist(a1, a2), math.dist(b1, b2)
        if not (0.5 < la < 0.015 * W and abs(la - lb) < 0.25 * max(la, lb)):
            continue
        flat = lambda p, q: abs(p.y - q.y) < 0.1 * math.dist(p, q)
        upright = lambda p, q: abs(p.x - q.x) < 0.1 * math.dist(p, q)
        if not ((flat(a1, a2) and upright(b1, b2)) or (upright(a1, a2) and flat(b1, b2))):
            continue
        ma, mb = ((a1.x + a2.x) / 2, (a1.y + a2.y) / 2), ((b1.x + b2.x) / 2, (b1.y + b2.y) / 2)
        if math.dist(ma, mb) <= 0.15 * la:
            found.append((ma[0] / W, ma[1] / H))
    return found


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
            r = page.rect
            labels = []
            for info in page.get_image_info():
                x1, y1, x2, y2 = info["bbox"]
                b = (x1 / r.width, y1 / r.height, x2 / r.width, y2 / r.height)
                if b[2] - b[0] < LABEL_IMAGE_MAX and b[3] - b[1] < LABEL_IMAGE_MAX:
                    labels.append(b)
            pages.append(RenderedPage(index=i, image=image, label_images=labels,
                                      plus_signs=_plus_signs(page)))

    return pages
