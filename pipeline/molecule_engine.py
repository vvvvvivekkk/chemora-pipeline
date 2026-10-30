"""
Stage 2a: image/PDF -> molecule SMILES  (OCSR = Optical Chemical Structure
Recognition — the chemistry version of OCR).

Two backends behind one idea:

* OpenChemIE (PyTorch) — the *preferred* path. It can take a WHOLE PDF and do its
  own figure-finding + recognition, which is more reliable than us slicing pages
  ourselves. So for this engine the orchestrator hands it the PDF directly
  (see `extract_molecules_from_pdf`). It also has an image path used as fallback.
* DECIMER (TensorFlow) — molecule-only. DECIMER recognises ONE molecule per
  image, so we segment a page into crops first, then recognise each crop
  (`recognize_molecules`).

API notes verified against the libraries' docs:
  - DECIMER:   predict_SMILES(image_input, confidence=False, hand_drawn=False)
               image_input may be a file path OR a numpy array; returns a SMILES
               string (or, with confidence=True, string + confidence info).
  - OpenChemIE: OpenChemIE(device=torch.device('cpu'|'cuda'))  # NOTE: torch.device
               .extract_molecules_from_figures_in_pdf(pdf_path)
               .extract_molecules_from_figures(list_of_images)   # fallback

Design rules (unchanged): lazy load, cache, never crash the app.
"""

from __future__ import annotations

import math
import re
import statistics
from dataclasses import dataclass

from PIL import Image

from . import config


@dataclass
class DetectedMolecule:
    """One molecule found in the document."""
    smiles: str
    confidence: float = 1.0
    page: int | None = None                 # 0-based page it was found on, if known
    bbox: tuple[int, int, int, int] | None = None
    crop: Image.Image | None = None
    molfile: str | None = None              # MolScribe's molblock: keeps 2D layout + R-labels


# Module-level cache so the (slow) model loads at most once per process.
_MODEL_CACHE: dict[str, object] = {}


# ---------------------------------------------------------------------------
# Shared OpenChemIE loader (used by this module AND reaction_engine.py)
# ---------------------------------------------------------------------------
def _get_openchemie():
    """Load (once) and return the OpenChemIE model, on GPU if available.

    IMPORTANT: OpenChemIE wants a torch.device, not the string 'cuda'/'cpu'.
    """
    if "openchemie" not in _MODEL_CACHE:
        import torch
        from openchemie import OpenChemIE  # type: ignore

        device = torch.device(config.resolve_device())
        _MODEL_CACHE["openchemie"] = OpenChemIE(device=device)
    return _MODEL_CACHE["openchemie"]


# ---------------------------------------------------------------------------
# Small helpers to normalise OpenChemIE's (nested, version-varying) output
# ---------------------------------------------------------------------------
def _coerce_smiles(entry) -> str | None:
    if isinstance(entry, str):
        return entry
    if isinstance(entry, dict):
        return entry.get("smiles") or entry.get("SMILES")
    return None


def _coerce_score(entry) -> float:
    # Prefer MolScribe's 0-1 recognition confidence; MolDetect's 'score' is not a
    # probability (values in the thousands), so it is only a last resort.
    if isinstance(entry, dict):
        try:
            return float(entry.get("confidence", entry.get("score", 1.0)) or 1.0)
        except (TypeError, ValueError):
            return 1.0
    return 1.0


def _molecules_from_openchemie_figure(fig: dict) -> list[DetectedMolecule]:
    """Pull molecules out of one OpenChemIE figure-result dict."""
    out: list[DetectedMolecule] = []
    # OpenChemIE reports 'page' as a 0-based index (pdfminer page_numbers), which
    # matches our RenderedPage.index directly — use it as-is (verified in source).
    page0 = fig.get("page") if isinstance(fig, dict) else None

    mol_items = []
    if isinstance(fig, dict):
        mol_items = fig.get("molecules") or fig.get("mols") or []
        # Some versions put a single molecule directly on the figure dict.
        if not mol_items and (fig.get("smiles") or fig.get("SMILES")):
            mol_items = [fig]

    for m in mol_items:
        smi = clean_model_smiles(_coerce_smiles(m))
        if smi:
            crop = m.get("image") if isinstance(m, dict) else None
            out.append(DetectedMolecule(
                smiles=smi,
                confidence=_coerce_score(m),
                page=page0,
                bbox=tuple(m["bbox"]) if isinstance(m, dict) and m.get("bbox") else None,
                # OpenChemIE crops are BGR numpy arrays; flip to RGB for display.
                crop=Image.fromarray(crop[:, :, ::-1]) if crop is not None and getattr(crop, "size", 0) else None,
                molfile=m.get("molfile") if isinstance(m, dict) else None,
            ))
    return out


def _pad_bbox(bbox, pad: float):
    """Grow a normalised (x1, y1, x2, y2) box by `pad` of its size, clipped to the page."""
    x1, y1, x2, y2 = bbox
    dx, dy = (x2 - x1) * pad, (y2 - y1) * pad
    return (max(0.0, x1 - dx), max(0.0, y1 - dy), min(1.0, x2 + dx), min(1.0, y2 + dy))


def postprocess_smiles(raw: dict, page=None) -> str | None:
    """All post-processing of one raw MolScribe result -> the SMILES shown on the canvas.

    page: the RenderedPage it came from (optional). Its picture-labels and drawn
    '+' signs are evidence for the fixes below.
    """
    mol, on_page = _molfile_on_page(raw)
    rgroups = _atoms_on_label_images(mol, on_page, page.label_images if page else ())
    cations = _atoms_at_plus_signs(mol, on_page, page.plus_signs if page else ())
    smiles = expand_unexpanded_labels(raw.get("smiles"), raw.get("molfile"), rgroups, cations)
    return clean_model_smiles(repair_bond_orders(smiles))


def _molfile_on_page(raw: dict):
    """MolScribe's molfile as a mol, plus a function: atom index -> (x, y) on the page."""
    molfile, bbox = raw.get("molfile"), raw.get("bbox")
    # crop image size; reaction parts carry only "size" (their crop isn't kept)
    size = raw["crop"].size if raw.get("crop") is not None else raw.get("size")
    if not (molfile and bbox and size and size[1]):
        return None, None
    from rdkit import Chem
    mol = Chem.MolFromMolBlock(molfile, sanitize=False, removeHs=False)
    if mol is None or not mol.GetNumConformers():
        return None, None
    conf, ratio = mol.GetConformer(), size[0] / size[1]
    x1, y1, x2, y2 = bbox

    def on_page(i):  # MolScribe writes (10*ratio*x, 1 - 10*y) for crop fractions x, y
        p = conf.GetAtomPosition(i)
        return x1 + p.x / (10 * ratio) * (x2 - x1), y1 + (1 - p.y) / 10 * (y2 - y1)

    return mol, on_page


def _atoms_at_plus_signs(mol, on_page, signs) -> set[int]:
    """Carbons drawn with a '+' next to them (carbocations); MolScribe ignores the sign.

    A sign is given to the nearest atom within one bond length, and only if that
    atom is a neutral carbon with three bonds' worth of bonding, so the '+'
    replaces its implicit H. A '+' over a reaction arrow (H+) is too far from
    any atom; an N next to a '+' is left alone (MolScribe reads N+ itself).
    """
    if not signs or mol is None or not mol.GetNumBonds():
        return set()
    bond = statistics.median(math.dist(on_page(b.GetBeginAtomIdx()), on_page(b.GetEndAtomIdx()))
                             for b in mol.GetBonds())
    found = set()
    for sign in signs:
        near = min(mol.GetAtoms(), key=lambda a: math.dist(on_page(a.GetIdx()), sign))
        if math.dist(on_page(near.GetIdx()), sign) > bond:
            continue
        if (near.GetSymbol() == "C" and near.GetFormalCharge() == 0
                and sum(b.GetBondTypeAsDouble() for b in near.GetBonds()) == 3):
            found.add(near.GetIdx())
    return found


def repair_bond_orders(smiles: str | None) -> str | None:
    """Fix a small double-bond mistake in an otherwise correctly read ring system.

    MolScribe can get every atom and bond of a fused aromatic system right but
    put two double bonds on one carbon (valence 5), so RDKit rejects the whole
    structure. If at most two atoms are over-valent, each by exactly one bond,
    the double bonds of the fused ring system around them are re-derived (marked
    aromatic, RDKit finds a valid alternating form). Atoms and connectivity are
    never changed. Anything more broken stays invalid, so the UI keeps warning
    the chemist instead of showing a valid-looking guess.
    """
    if not smiles:
        return smiles
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    if Chem.MolFromSmiles(smiles) is not None:
        return smiles
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return smiles
    mol.UpdatePropertyCache(strict=False)
    table = Chem.GetPeriodicTable()
    over = []
    for a in mol.GetAtoms():
        if a.GetFormalCharge() or a.GetAtomicNum() == 0:
            continue
        excess = sum(b.GetBondTypeAsDouble() for b in a.GetBonds()) - max(table.GetValenceList(a.GetAtomicNum()))
        if excess > 0:
            over.append((a.GetIdx(), excess))
    if not over or len(over) > 2 or any(e != 1 for _, e in over):
        return smiles
    Chem.FastFindRings(mol)
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    system = {i for idx, _ in over for r in rings if idx in r for i in r}
    grew = True
    while grew:  # the whole fused ring system, not just the first ring
        grew = False
        for r in rings:
            if r & system and not r <= system:
                system |= r
                grew = True
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in system and j in system and any(i in r and j in r for r in rings):
            b.SetBondType(Chem.BondType.AROMATIC)
            b.SetIsAromatic(True)
    for i in system:
        mol.GetAtomWithIdx(i).SetIsAromatic(True)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return smiles
    return Chem.MolToSmiles(mol)


# What a drawn structure can plainly end in. Anything else at a chain end that
# sits on a label picture (Pd, Eu, Rf, deuterium...) is MolScribe guessing.
_PLAIN_ENDS = {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B", "Si", "Se", "H"}


def _atoms_on_label_images(mol, on_page, boxes) -> set[int]:
    """Atoms MolScribe placed on a picture-label (an R', R'' drawn as an image).

    MolScribe can't read such a picture and emits an arbitrary element there.
    For each glyph-sized image near the crop we take the nearest chain-end atom
    and mark it as an R-group only if it is not something a structure plainly
    ends in: a real CH3, F or OH next to an image is never touched.
    """
    if not boxes or mol is None:
        return set()
    ends = [a for a in mol.GetAtoms() if a.GetDegree() == 1]
    found = set()
    for bx1, by1, bx2, by2 in boxes:
        centre, reach = ((bx1 + bx2) / 2, (by1 + by2) / 2), 1.5 * max(bx2 - bx1, by2 - by1)
        near = min(ends, key=lambda a: math.dist(on_page(a.GetIdx()), centre), default=None)
        if near is None or math.dist(on_page(near.GetIdx()), centre) > reach:
            continue
        if near.GetAtomicNum() and (near.GetSymbol() not in _PLAIN_ENDS or near.GetIsotope()):
            found.add(near.GetIdx())
    return found


# 'NMe2', 'NEt2', ...: an element followed by a group and a count
_GROUP_TIMES_N = re.compile(r"^([A-Z][a-z]?)(Me|Et|Pr|Bu|Ph)(\d)$")


def expand_unexpanded_labels(smiles: str | None, molfile: str | None,
                             rgroup_atoms=frozenset(), cation_atoms=frozenset()) -> str | None:
    """Second chance for atom labels MolScribe read but could not expand.

    MolScribe left those atoms as '*'. Its formula parser reads 'NMe2' greedily
    as 'NMe' x 2 ('NMe' is one of its own abbreviations) and gives up, although
    it parses 'N(Me)2' fine. We redo MolScribe's own expansion from its molfile
    (which keeps the labels) with that retry added, and use the result only if
    it has fewer '*' atoms than before — otherwise the original stands.
    rgroup_atoms / cation_atoms (molfile atom indices) are first turned into
    R-groups / given a + charge; then the rebuilt result is used as is.
    """
    if not smiles or not molfile or ("*" not in smiles and not rgroup_atoms and not cation_atoms):
        return smiles
    try:
        from rdkit import Chem
        from molscribe import chemistry  # type: ignore
    except Exception:
        return smiles
    mol = Chem.MolFromMolBlock(molfile, sanitize=False, removeHs=False)
    if mol is None:
        return smiles
    for i in rgroup_atoms:
        atom = mol.GetAtomWithIdx(i)
        atom.SetAtomicNum(0)
        atom.SetIsotope(0)
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(True)
        Chem.SetAtomAlias(atom, "R")  # MolScribe leaves R-groups as '*'
    for i in cation_atoms:
        atom = mol.GetAtomWithIdx(i)
        atom.SetFormalCharge(1)  # the '+' takes the place of the implicit H
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            # RDKit names molfile dummies 'R'; MolScribe only expands atoms named '*'
            if atom.HasProp("dummyLabel"):
                atom.ClearProp("dummyLabel")
            if atom.HasProp("molFileAlias"):
                Chem.SetAtomAlias(atom, atom.GetProp("molFileAlias"))

    parse = chemistry.get_smiles_from_symbol

    def parse_with_retry(symbol, mol, atom, bonds):
        out = parse(symbol, mol, atom, bonds)
        m = _GROUP_TIMES_N.match(symbol) if out is None and isinstance(symbol, str) else None
        return parse(f"{m[1]}({m[2]}){m[3]}", mol, atom, bonds) if m else out

    chemistry.get_smiles_from_symbol = parse_with_retry  # only for this call
    try:
        redone, _ = chemistry._expand_functional_group(mol, {})
    except Exception:
        return smiles
    finally:
        chemistry.get_smiles_from_symbol = parse
    if rgroup_atoms or cation_atoms:
        return redone
    return redone if redone.count("*") < smiles.count("*") else smiles


def clean_model_smiles(smiles: str | None) -> str | None:
    """Undo systematic MolScribe misreads seen on real papers; None = drop it.

    * R-group labels drawn in boxes/circles (R', R'') come back as heavy,
      radioactive elements (Rn, Pa, Rf, Fl, ... seen on real papers) that never
      occur in an organic chemistry paper -> R-group *.
    * A crop holds one molecule, so extra neutral fragments are leakage from
      nearby text/catalyst labels (Pd, W, CH4), logos (Co.Co.Co) or the edge of
      a neighbouring structure -> only the largest fragment is kept, plus any
      charged ions (Na+, Cl-) so real salts survive.
    * Nothing carbon-based left (a logo read as '*.*') -> not a structure.
    Unparseable SMILES are returned unchanged so the UI still flags them.
    """
    if not smiles:
        return None
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
    except Exception:
        return smiles
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles
    rw = Chem.RWMol(mol)
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() >= 84:           # Po, At, Rn, ... Pa, U, ... Rf, Fl
            atom.SetAtomicNum(0)
            atom.SetIsotope(0)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)
            atom.SetNumRadicalElectrons(0)
            atom.SetFormalCharge(0)
    frags = sorted(Chem.GetMolFrags(rw.GetMol(), asMols=True, sanitizeFrags=False),
                   key=lambda f: f.GetNumHeavyAtoms(), reverse=True)
    keep = [f for i, f in enumerate(frags)
            if (i == 0 and f.GetNumHeavyAtoms() > 1)
            or any(a.GetFormalCharge() for a in f.GetAtoms())]
    if not any(a.GetAtomicNum() == 6 for f in keep for a in f.GetAtoms()):
        return None
    out = keep[0]
    for f in keep[1:]:
        out = Chem.CombineMols(out, f)
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return smiles
    return Chem.MolToSmiles(out)


# ---------------------------------------------------------------------------
# PREFERRED PATH: whole-PDF extraction with OpenChemIE
# ---------------------------------------------------------------------------
def extract_molecules_from_pdf(pdf_path: str) -> list[DetectedMolecule]:
    """Run OpenChemIE over a whole PDF and return every molecule found.

    Each DetectedMolecule carries the page it was found on (when OpenChemIE
    reports it) so the orchestrator can group results per page.
    Returns [] on any failure — never raises.
    """
    num_pages = None if config.MAX_PAGES == 0 else config.MAX_PAGES
    try:
        model = _get_openchemie()
        results = model.extract_molecules_from_figures_in_pdf(pdf_path, num_pages=num_pages)
    except Exception as exc:
        return [DetectedMolecule(smiles=f"[molecule model error: {exc}]", confidence=0.0)]

    molecules: list[DetectedMolecule] = []
    for fig in results or []:
        molecules.extend(_molecules_from_openchemie_figure(fig))
    return molecules


# ---------------------------------------------------------------------------
# PAGE-IMAGE PATH: molecules from pre-rendered page images (no layoutparser)
# ---------------------------------------------------------------------------
def extract_molecules_from_images(images: list[Image.Image], pages=None) -> list[DetectedMolecule]:
    """Run OpenChemIE's MolDetect+MolScribe directly on page images.

    We pass the pages we already rendered (PyMuPDF) instead of the PDF, so
    OpenChemIE never calls ``extract_figures_from_pdf`` -> layoutparser, whose
    PubLayNet weights are hosted on a now-broken Dropbox link (returns HTML).
    MolDetect finds molecule bounding boxes within each page image on its own.

    Results come back aligned to the input list, so the image's position is its
    page index. Returns [] / an error marker on failure — never raises.
    """
    try:
        raw = read_molecules_raw(images)
    except Exception as exc:
        return [DetectedMolecule(smiles=f"[molecule model error: {exc}]", confidence=0.0)]
    return molecules_from_raw(raw, pages)


def read_molecules_raw(images: list[Image.Image]) -> list[list[dict]]:
    """Model step only (slow). Per page: [{bbox, smiles, molfile, confidence,
    crop}] exactly as MolDetect + MolScribe produced them, before any cleanup.
    bench/ caches this so post-processing changes can be re-scored without
    re-running the model. Raises on model failure.
    """
    from openchemie.utils import clean_bbox_output, convert_to_cv2  # type: ignore
    model = _get_openchemie()
    # Same steps as OpenChemIE.extract_molecules_from_figures, unrolled because
    # MolScribe crashes on an empty batch ("not enough values to unpack"):
    # a document with no structures must give [] rather than an error.
    bboxes = model.extract_molecule_bboxes_from_figures(images, batch_size=4)
    for page_boxes in bboxes:
        for b in page_boxes:
            b["bbox"] = _pad_bbox(b["bbox"], config.CROP_PAD)
    figures = [convert_to_cv2(im) for im in images]
    results, crops, refs = clean_bbox_output(figures, bboxes)
    if crops:
        infos = model.molscribe.predict_images(crops, return_confidence=True, batch_size=4)
        for info, ref in zip(infos, refs):
            ref.update(info)
    pages = []
    for fig in results:
        items = []
        for m in fig.get("molecules", []):
            crop = m.get("image")
            items.append({
                "bbox": [float(v) for v in m["bbox"]],
                "smiles": _coerce_smiles(m),
                "molfile": m.get("molfile"),
                "confidence": _coerce_score(m),
                # OpenChemIE crops are BGR numpy arrays; flip to RGB for display.
                "crop": Image.fromarray(crop[:, :, ::-1]) if crop is not None and getattr(crop, "size", 0) else None,
            })
        pages.append(items)
    return pages


def molecules_from_raw(raw_pages: list[list[dict]], pages=None) -> list[DetectedMolecule]:
    """Post-processing step (fast): cleanup -> DetectedMolecule. Page = list index.

    pages: the RenderedPages (optional), for evidence from the PDF itself.
    """
    molecules: list[DetectedMolecule] = []
    for page_idx, items in enumerate(raw_pages):
        page = pages[page_idx] if pages and page_idx < len(pages) else None
        for m in items:
            smi = postprocess_smiles(m, page)
            if smi:
                molecules.append(DetectedMolecule(
                    smiles=smi, confidence=m.get("confidence", 1.0), page=page_idx,
                    bbox=tuple(m["bbox"]) if m.get("bbox") else None,
                    crop=m.get("crop"), molfile=m.get("molfile"),
                ))
    return molecules


# ---------------------------------------------------------------------------
# PER-IMAGE PATH: used by the DECIMER engine (and as an OpenChemIE fallback)
# ---------------------------------------------------------------------------
def recognize_molecules(page_image: Image.Image) -> list[DetectedMolecule]:
    """Find molecules on a single page image. Used for the DECIMER engine.

    The orchestrator only calls this for the decimer engine.
    """
    try:
        return _recognize_with_decimer(page_image)
    except Exception as exc:
        return [DetectedMolecule(smiles=f"[recognition error: {exc}]", confidence=0.0)]


def _recognize_with_decimer(page_image: Image.Image) -> list[DetectedMolecule]:
    """DECIMER path: segment the page into structure crops, recognise each."""
    import numpy as np

    # Step 1: segment the page into individual structure images.
    crops: list[Image.Image] = []
    try:
        from decimer_segmentation import segment_chemical_structures  # type: ignore
        segments = segment_chemical_structures(np.array(page_image))
        crops = [Image.fromarray(s) for s in segments]
    except Exception:
        crops = [page_image]  # no segmenter -> treat whole page as one structure

    # Step 2: recognise each crop. DECIMER accepts a numpy array directly, so no
    # temp files are needed. confidence=True gives us a real confidence value.
    from DECIMER import predict_SMILES  # type: ignore

    molecules: list[DetectedMolecule] = []
    for crop in crops:
        try:
            out = predict_SMILES(np.array(crop), confidence=True)
            smi, conf = _split_decimer_output(out)
            if smi:
                molecules.append(DetectedMolecule(smiles=smi, confidence=conf, crop=crop))
        except Exception:
            continue
    return molecules


def _split_decimer_output(out) -> tuple[str, float]:
    """DECIMER with confidence=True returns either a (smiles, conf) pair, a list
    of (token, conf) tuples, or just a string depending on version. Normalise."""
    if isinstance(out, str):
        return out, 1.0
    if isinstance(out, tuple) and len(out) == 2 and isinstance(out[0], str):
        try:
            return out[0], float(out[1])
        except (TypeError, ValueError):
            return out[0], 1.0
    if isinstance(out, (list, tuple)):
        # list of (token, confidence) -> join tokens, average confidence
        tokens, confs = [], []
        for item in out:
            if isinstance(item, (list, tuple)) and len(item) == 2:
                tokens.append(str(item[0]))
                try:
                    confs.append(float(item[1]))
                except (TypeError, ValueError):
                    pass
        smi = "".join(tokens)
        conf = sum(confs) / len(confs) if confs else 1.0
        return smi, conf
    return str(out), 1.0
