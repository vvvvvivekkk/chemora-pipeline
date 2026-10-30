"""
The orchestrator: the ONE function the UI calls (`process_pdf`).

It wires the stages together and chooses the right strategy per engine:

  * openchemie : hand the WHOLE PDF to OpenChemIE (it finds figures itself),
                 then map the molecules/reactions back onto the rendered pages by
                 page number. This is the most reliable path.
  * decimer    : go page-by-page — segment each page image, recognise each crop.
  * manual     : just render the pages; the user types SMILES in the UI.

Everything is wrapped so recognition problems are recorded in `.status`, never
raised — the UI always has something to show.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from PIL import Image

from . import config, pdf_utils
from .molecule_engine import (
    recognize_molecules,
    extract_molecules_from_images,
    DetectedMolecule,
)
from .reaction_engine import extract_reactions_from_images, DetectedReaction
from .validate import validate_smiles, ValidatedMolecule


@dataclass
class PageResult:
    index: int
    image: Image.Image
    molecules: list[ValidatedMolecule] = field(default_factory=list)
    reactions: list[DetectedReaction] = field(default_factory=list)


@dataclass
class ExtractionResult:
    pages: list[PageResult] = field(default_factory=list)
    status: list[str] = field(default_factory=list)
    engine: str = config.ENGINE

    @property
    def total_molecules(self) -> int:
        return sum(len(p.molecules) for p in self.pages)

    @property
    def total_reactions(self) -> int:
        return sum(len(p.reactions) for p in self.pages)


def process_pdf(pdf_bytes: bytes, progress_callback=None, engine: str | None = None) -> ExtractionResult:
    # The engine is per call, never written back to `config`: config is shared by
    # every browser session on the server, so mutating it leaks one user's choice.
    engine = engine or config.ENGINE
    result = ExtractionResult(engine=engine)

    def _progress(frac, msg):
        result.status.append(msg)
        if progress_callback:
            progress_callback(frac, msg)

    # --- Stage 1: render pages (needed for display in every mode) ---------------
    try:
        pages = pdf_utils.pdf_to_images(pdf_bytes)
    except Exception as exc:
        result.status.append(f"Could not read PDF: {exc}")
        return result
    if not pages:
        result.status.append("No pages found in PDF.")
        return result

    for p in pages:
        result.pages.append(PageResult(index=p.index, image=p.image))
    by_index = {pr.index: pr for pr in result.pages}
    _progress(0.1, f"Rendered {len(pages)} page(s) at {config.PDF_RENDER_DPI} DPI.")

    # --- MANUAL: stop here ------------------------------------------------------
    if engine == "manual":
        _progress(1.0, "Manual mode: no model run. Use the canvas to enter SMILES.")
        return result

    # --- OPENCHEMIE: recognise from the rendered page images -------------------
    # We feed the page images we already rendered (above) straight to OpenChemIE's
    # MolDetect/MolScribe/RxnScribe. This deliberately avoids OpenChemIE's own
    # PDF path (extract_*_from_figures_in_pdf), which runs layoutparser to find
    # figures first — and layoutparser's PubLayNet weights are hosted on a broken
    # Dropbox link that returns HTML on Windows. MolDetect locates structures in
    # each full page image on its own, so figure-finding isn't needed.
    if engine == "openchemie":
        images = [pr.image for pr in result.pages]

        _progress(0.2, "Loading models & detecting molecules in page images…")
        for det in extract_molecules_from_images(images, [p.label_images for p in pages]):
            _attach_molecule(by_index, result.pages, det)

        _progress(0.6, "Detecting reactions in page images…")
        for rxn in extract_reactions_from_images(images, [p.label_images for p in pages]):
            _attach_reaction(by_index, result.pages, rxn)

        _progress(
            1.0,
            f"Done. {result.total_molecules} molecule(s) and "
            f"{result.total_reactions} reaction(s).",
        )
        return result

    # --- DECIMER: page-by-page --------------------------------------------------
    n = len(pages)
    for i, p in enumerate(pages):
        pr = by_index[p.index]
        for det in recognize_molecules(p.image):
            pr.molecules.append(_validate_detection(det))
        frac = 0.1 + 0.9 * (i + 1) / n
        _progress(frac, f"Page {p.index + 1}/{n}: {len(pr.molecules)} molecule(s).")

    _progress(1.0, f"Done. {result.total_molecules} molecule(s) across {n} page(s).")
    return result


# ---------------------------------------------------------------------------
# helpers to place a detection on the right page and validate it
# ---------------------------------------------------------------------------
def _validate_detection(det: DetectedMolecule) -> ValidatedMolecule:
    v = validate_smiles(det.smiles)
    if det.confidence < config.LOW_CONFIDENCE_THRESHOLD and v.is_valid:
        v.message += f"  (low model confidence: {det.confidence:.2f})"
    return v


def _attach_molecule(by_index, pages, det: DetectedMolecule):
    """Place a molecule on its reported page, falling back to the first page."""
    target = by_index.get(det.page) if det.page is not None else None
    if target is None:
        target = pages[0]
    target.molecules.append(_validate_detection(det))


def _attach_reaction(by_index, pages, rxn: DetectedReaction):
    target = by_index.get(rxn.page) if rxn.page is not None else None
    if target is None:
        target = pages[0]
    target.reactions.append(rxn)
