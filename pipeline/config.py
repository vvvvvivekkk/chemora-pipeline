"""
Central configuration for the Chemora Extract pipeline.

Everything tunable lives here so you never have to hunt through the code to
change behaviour. The UI and every pipeline stage read their settings from this
one file.
"""

import os
import sys


def _ensure_poppler_on_path() -> None:
    """Make the conda env's poppler binaries findable (Windows).

    OpenChemIE uses pdf2image, which shells out to poppler's `pdfinfo`/`pdftoppm`.
    conda-forge installs those into `<env>\\Library\\bin`, which is only on PATH
    when the env is *activated*. If the app is launched with the interpreter
    directly (e.g. a VS Code run config, or `python app.py` without activation),
    PATH is missing that dir and OpenChemIE fails with "Unable to get page count.
    Is poppler installed and in PATH?". Prepending it here makes extraction work
    no matter how the app was launched. No-op off Windows / if the dir is absent.
    """
    if os.name != "nt":
        return
    lib_bin = os.path.join(sys.prefix, "Library", "bin")
    if os.path.isdir(lib_bin) and lib_bin.lower() not in os.environ.get("PATH", "").lower():
        os.environ["PATH"] = lib_bin + os.pathsep + os.environ.get("PATH", "")


_ensure_poppler_on_path()

# Windows compatibility shims for the OpenChemIE dependency stack (no-op elsewhere).
try:
    from . import _win_compat
    _win_compat.apply()
except Exception:
    pass


# ---------------------------------------------------------------------------
# WHICH RECOGNITION ENGINE TO USE
# ---------------------------------------------------------------------------
# "openchemie" -> molecules AND reactions (PyTorch).  Closest to the full vision.
# "decimer"    -> molecules only (TensorFlow). Simpler fallback.
# "manual"     -> no model loads at all; you paste SMILES and see it drawn.
#                 Always works, even with no GPU. Great for demoing the UI.
#
# You can also override this without editing the file by setting an environment
# variable before launching, e.g.:   export CHEMORA_ENGINE=decimer
# ---------------------------------------------------------------------------
ENGINE = os.environ.get("CHEMORA_ENGINE", "openchemie")


# ---------------------------------------------------------------------------
# PDF RENDERING
# ---------------------------------------------------------------------------
# How sharp to render each PDF page before the models look at it.
# Higher DPI = clearer structures = better recognition, but bigger images and
# slower processing. 300 is a good balance for printed journal figures.
PDF_RENDER_DPI = int(os.environ.get("CHEMORA_PDF_DPI", "300"))

# Safety cap so a 200-page PDF can't hang the demo. Only the first N pages are
# processed. Set to 0 for "all pages".
MAX_PAGES = int(os.environ.get("CHEMORA_MAX_PAGES", "10"))


# ---------------------------------------------------------------------------
# DEVICE
# ---------------------------------------------------------------------------
# "cuda" uses a GPU if one is present; "cpu" forces CPU. "auto" picks the best
# available at runtime (see _resolve_device below).
DEVICE = os.environ.get("CHEMORA_DEVICE", "auto")


def resolve_device() -> str:
    """Return 'cuda' or 'cpu', honouring the DEVICE setting.

    We import torch lazily inside the function so that the 'manual' engine, which
    needs no deep-learning framework at all, never pays the cost of importing it.
    """
    if DEVICE in ("cuda", "cpu"):
        return DEVICE
    try:
        import torch
        return "cuda" if torch.cuda.is_available() else "cpu"
    except Exception:
        # torch not installed (e.g. manual mode) -> CPU is the only honest answer.
        return "cpu"


# ---------------------------------------------------------------------------
# RECOGNITION CONFIDENCE
# ---------------------------------------------------------------------------
# Results below this confidence get flagged in the UI as "please double-check".
# This is a UI hint, not a hard filter — we never silently drop a result, because
# a chemist deciding what to trust is the whole point of human-in-the-loop.
LOW_CONFIDENCE_THRESHOLD = float(os.environ.get("CHEMORA_LOW_CONF", "0.5"))


# ---------------------------------------------------------------------------
# CROP PADDING
# ---------------------------------------------------------------------------
# Grow each MolDetect box by this fraction of its width/height before MolScribe
# reads it. Off by default: on real papers (bench/benchmark.py papers) 0.08 did
# not recover labels at the box edge (NMe2 still read as R) and pulled captions
# and neighbouring structures into crops. Kept as a knob for other layouts.
CROP_PAD = float(os.environ.get("CHEMORA_CROP_PAD", "0"))
