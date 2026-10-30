"""
Stage 3 of the pipeline: validate and draw a SMILES with RDKit.

A recognition model can output a SMILES string that isn't real chemistry (wrong
valence, nonsense symbols, a typo in the ring). RDKit is the standard open-source
"chemistry engine" — we use it to:

  1. CHECK the SMILES actually parses into a valid molecule.
  2. CANONICALIZE it (rewrite it in one standard form, so "the same molecule"
     always looks the same as text).
  3. DRAW a clean picture of it (so the UI can show the result even before the
     interactive canvas loads).

SMILES = "Simplified Molecular-Input Line-Entry System": a molecule written as a
short line of text. Example: aspirin is  CC(=O)Oc1ccccc1C(=O)O
"""

from __future__ import annotations

import io
from dataclasses import dataclass

from PIL import Image


@dataclass
class ValidatedMolecule:
    """The result of checking one SMILES string."""
    raw_smiles: str              # exactly what the model produced
    canonical_smiles: str | None # RDKit's standard form, or None if invalid
    is_valid: bool               # did it parse as real chemistry?
    message: str                 # human-readable note (e.g. why it was rejected)
    image: Image.Image | None    # a drawn picture of the molecule, or None


def validate_smiles(smiles: str) -> ValidatedMolecule:
    """Validate a single SMILES string and render it to an image."""
    smiles = (smiles or "").strip()

    if not smiles:
        return ValidatedMolecule(smiles, None, False, "Empty SMILES.", None)

    # Imported lazily: RDKit is a big dependency, and keeping the import inside the
    # function means a caller that only needs other utilities doesn't pay for it.
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
        from rdkit.Chem.Draw import rdMolDraw2D  # noqa: F401 (import checks availability)
    except Exception as exc:  # RDKit missing
        return ValidatedMolecule(
            smiles, None, False,
            f"RDKit not available ({exc}). Cannot validate.", None,
        )

    # Chem.MolFromSmiles returns None if the string isn't a valid molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ValidatedMolecule(
            smiles, None, False,
            "RDKit could not parse this as a valid molecule — likely a "
            "recognition error. Edit it on the canvas.",
            None,
        )

    canonical = Chem.MolToSmiles(mol)  # the one standard spelling

    # Draw a clean 2D depiction.
    try:
        pil_img = Draw.MolToImage(mol, size=(350, 300))
    except Exception:
        pil_img = None

    return ValidatedMolecule(
        raw_smiles=smiles,
        canonical_smiles=canonical,
        is_valid=True,
        message="Valid structure.",
        image=pil_img,
    )


def smiles_to_png_bytes(smiles: str, size: tuple[int, int] = (350, 300)) -> bytes | None:
    """Convenience: return a PNG (as bytes) for a SMILES, or None if it fails."""
    result = validate_smiles(smiles)
    if result.image is None:
        return None
    buf = io.BytesIO()
    result.image.save(buf, format="PNG")
    return buf.getvalue()
