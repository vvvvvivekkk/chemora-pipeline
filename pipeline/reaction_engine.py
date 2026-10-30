"""
Stage 2b: reaction diagram -> structured reaction (reactants / conditions /
products), using OpenChemIE (which bundles RxnScribe for this).

Verified API:
  model.extract_reactions_from_figures_in_pdf(pdf_path)   # preferred (whole PDF)
  model.extract_reactions_from_figures(list_of_images)    # fallback (per image)

Only the 'openchemie' engine supports reactions; 'decimer' and 'manual' return [].
Same rules: lazy load, cache, never crash.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from PIL import Image

from . import config
from .molecule_engine import _get_openchemie, _coerce_smiles


@dataclass
class DetectedReaction:
    """One reaction scheme found in the document."""
    reactants: list[str] = field(default_factory=list)   # SMILES
    products: list[str] = field(default_factory=list)     # SMILES
    conditions: list[str] = field(default_factory=list)   # free text
    confidence: float = 1.0
    page: int | None = None                               # 0-based page, if known

    def as_reaction_smiles(self) -> str:
        """Render as reaction SMILES: 'reactant.reactant>>product.product'."""
        return f"{'.'.join(self.reactants)}>>{'.'.join(self.products)}"


def extract_reactions_from_pdf(pdf_path: str) -> list[DetectedReaction]:
    """Run OpenChemIE reaction extraction over a whole PDF. Returns [] on failure."""
    if config.ENGINE != "openchemie":
        return []
    num_pages = None if config.MAX_PAGES == 0 else config.MAX_PAGES
    try:
        model = _get_openchemie()
        results = model.extract_reactions_from_figures_in_pdf(pdf_path, num_pages=num_pages)
    except Exception as exc:
        return [DetectedReaction(conditions=[f"[reaction model error: {exc}]"],
                                 confidence=0.0)]

    reactions: list[DetectedReaction] = []
    for fig in results or []:
        reactions.extend(_reactions_from_openchemie_figure(fig))
    return reactions


def extract_reactions_from_images(images: list[Image.Image]) -> list[DetectedReaction]:
    """Run OpenChemIE's RxnScribe directly on pre-rendered page images.

    Same rationale as ``molecule_engine.extract_molecules_from_images``: bypass
    the PDF -> layoutparser figure finder (broken model download on Windows) by
    feeding the page images we already have. Results are aligned to the input
    list, so an image's position is its page index. The orchestrator only calls
    this for the openchemie engine.
    """
    try:
        model = _get_openchemie()
        results = model.extract_reactions_from_figures(images)
    except Exception as exc:
        return [DetectedReaction(conditions=[f"[reaction model error: {exc}]"],
                                 confidence=0.0)]

    reactions: list[DetectedReaction] = []
    for page_idx, fig in enumerate(results or []):
        for rxn in _reactions_from_openchemie_figure(fig):
            rxn.page = page_idx
            reactions.append(rxn)
    return reactions


def _reactions_from_openchemie_figure(fig: dict) -> list[DetectedReaction]:
    """Pull reactions out of one OpenChemIE figure-result dict."""
    out: list[DetectedReaction] = []
    # OpenChemIE 'page' is 0-based (verified in source) — matches our page index.
    page0 = fig.get("page") if isinstance(fig, dict) else None

    rxn_items = fig.get("reactions", []) if isinstance(fig, dict) else []
    # Some versions put a single reaction's keys directly on the figure dict.
    if not rxn_items and isinstance(fig, dict) and (
        "reactants" in fig or "products" in fig
    ):
        rxn_items = [fig]

    for rxn in rxn_items:
        if not isinstance(rxn, dict):
            continue
        out.append(DetectedReaction(
            reactants=_smiles_list(rxn.get("reactants")),
            products=_smiles_list(rxn.get("products")),
            conditions=_text_list(rxn.get("conditions")),
            confidence=_score(rxn),
            page=page0,
        ))
    return out


def _smiles_list(items) -> list[str]:
    out = []
    for it in items or []:
        smi = _coerce_smiles(it)
        if smi:
            out.append(smi)
    return out


def _text_list(items) -> list[str]:
    """Flatten condition entries to strings.

    OpenChemIE conditions look like {'category':..., 'text': [list of str]} — note
    'text' is a LIST (verified in source), so we extend rather than append. We also
    accept a few other shapes defensively.
    """
    out: list[str] = []
    for it in items or []:
        if isinstance(it, str):
            out.append(it)
        elif isinstance(it, dict):
            txt = it.get("text")
            if isinstance(txt, list):
                out.extend(str(t) for t in txt if t)
            elif isinstance(txt, str) and txt:
                out.append(txt)
            else:
                alt = it.get("label") or it.get("smiles") or it.get("SMILES")
                if alt:
                    out.append(alt)
    return out


def _score(rxn: dict) -> float:
    try:
        return float(rxn.get("score", rxn.get("confidence", 1.0)) or 1.0)
    except (TypeError, ValueError):
        return 1.0
