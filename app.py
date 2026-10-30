"""
Chemora Extract — Streamlit UI.

The whole app is one screen with three parts:
  1. Sidebar  : upload a PDF, pick the engine, run.
  2. Left     : the original PDF page (so the chemist sees the source).
  3. Right    : the extracted structures on an EDITABLE Ketcher canvas
                (so the chemist reviews and fixes them).

The UI never touches a model directly — it only calls
`pipeline.orchestrator.process_pdf`. That keeps this file about *presentation*,
and keeps the chemistry/ML in the pipeline package where it can be reused.
"""

from __future__ import annotations

import csv
import io

import streamlit as st

from pipeline import config
from pipeline.orchestrator import process_pdf, ExtractionResult
from pipeline.validate import validate_smiles

# The interactive chemical canvas. Imported defensively so the app still loads
# (and can explain itself) even if the component isn't installed yet.
try:
    from streamlit_ketcher import st_ketcher
    _KETCHER_OK = True
except Exception:
    _KETCHER_OK = False


# ---------------------------------------------------------------------------
# Page setup
# ---------------------------------------------------------------------------
st.set_page_config(page_title="Chemora Extract", page_icon="🧪", layout="wide")

st.title("🧪 Chemora Extract")
st.caption(
    "Drop a chemistry PDF → extract its structures and reactions → "
    "review and correct them on an editable canvas."
)


# ---------------------------------------------------------------------------
# Sidebar: controls
# ---------------------------------------------------------------------------
with st.sidebar:
    st.header("1 · Upload")
    uploaded = st.file_uploader("Chemistry PDF", type=["pdf"])

    st.header("2 · Engine")
    engine_choice = st.selectbox(
        "Recognition engine",
        options=["openchemie", "decimer", "manual"],
        index=["openchemie", "decimer", "manual"].index(config.ENGINE)
        if config.ENGINE in ("openchemie", "decimer", "manual") else 0,
        help=(
            "openchemie = molecules + reactions · "
            "decimer = molecules only · "
            "manual = no model, you type SMILES (always works)."
        ),
    )
    run = st.button("Run extraction", type="primary", use_container_width=True)
    # Filled with the export button once results (and the chemist's edits) exist.
    export_slot = st.empty()

    st.divider()
    st.caption(
        "⚠️ First run downloads model weights and can be slow on CPU. "
        "For a live demo, run once beforehand so weights are cached, or use a "
        "GPU Space/Colab."
    )
    if not _KETCHER_OK:
        st.warning(
            "`streamlit-ketcher` isn't installed, so the editable canvas is "
            "disabled — structures will still show as static images. "
            "Add `streamlit-ketcher` to requirements.txt to enable editing."
        )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _render_canvas(smiles: str, key: str) -> str:
    """Show a SMILES on the editable Ketcher canvas (or a static image) and
    return the structure as it stands after the chemist's edits."""
    if _KETCHER_OK:
        # st_ketcher returns the seed until the chemist clicks Apply, then
        # whatever they drew — that's how the correction flows back.
        return st_ketcher(smiles, key=key) or ""
    # Fallback: static picture via RDKit (no editing possible).
    v = validate_smiles(smiles)
    if v.image is not None:
        st.image(v.image, caption=smiles)
    else:
        st.code(smiles)
    return smiles


def _canonical(smiles: str) -> str:
    return validate_smiles(smiles).canonical_smiles or smiles


def _review(seed: str, key: str, page_no: int, source: str, rows: list[dict]):
    """Canvas for one structure; records the chemist's final version in `rows`."""
    final = _render_canvas(seed, key)
    if final and seed and _canonical(final) == _canonical(seed):
        final = seed  # Ketcher re-spells SMILES; only a real change is an edit.
    if final != seed:
        if not final:
            st.caption("🗑️ Cleared — left out of the export.")
            return
        v = validate_smiles(final)
        if v.is_valid:
            st.success(f"✏️ Edited → `{v.canonical_smiles}`")
        else:
            st.error(f"✏️ Edited structure is not valid: {v.message}")
        final, source = v.canonical_smiles or final, f"{source}+edited"
    if final:
        rows.append({"page": page_no, "source": source, "smiles": final,
                     "valid": validate_smiles(final).is_valid})


def _export_csv(rows: list[dict]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=["page", "source", "smiles", "valid"])
    w.writeheader()
    w.writerows(rows)
    return buf.getvalue()


def _show_results(result: ExtractionResult):
    """Render the full result: a summary, then each page side-by-side."""
    # Top-line summary.
    c1, c2, c3 = st.columns(3)
    c1.metric("Pages", len(result.pages))
    c2.metric("Molecules", result.total_molecules)
    c3.metric("Reactions", result.total_reactions)

    with st.expander("Run log", expanded=False):
        for line in result.status:
            st.write("•", line)

    # The reviewed set: what the chemist ends up with after edits and additions.
    rows: list[dict] = []

    # Per-page: original on the left, extracted + canvas on the right.
    for page in result.pages:
        page_no = page.index + 1
        st.divider()
        st.subheader(f"Page {page_no}")
        left, right = st.columns([1, 1])

        with left:
            st.image(page.image, caption="Original page", use_container_width=True)

        with right:
            if result.engine == "manual":
                st.info("Manual mode — no model ran. Add structures below.")
            elif not page.molecules and not page.reactions:
                st.info("No structures or reactions detected on this page.")

            # Molecules
            for mi, mol in enumerate(page.molecules):
                label = mol.canonical_smiles or mol.raw_smiles
                badge = "✅" if mol.is_valid else "⚠️"
                st.markdown(f"**{badge} Molecule {mi + 1}** — `{label}`")
                if not mol.is_valid:
                    st.caption(mol.message)
                # Editable canvas, seeded with the recognised structure. An
                # unparseable seed would leave Ketcher blank and stuck, so the
                # chemist draws it from scratch instead.
                _review(mol.canonical_smiles or "", f"canvas_{page.index}_{mi}",
                        page_no, "extracted", rows)

            # Reactions
            for ri, rxn in enumerate(page.reactions):
                st.markdown(f"**⚗️ Reaction {ri + 1}**")
                st.code(rxn.as_reaction_smiles(), language="text")
                if rxn.conditions:
                    st.caption("Conditions: " + "; ".join(rxn.conditions))
                rows.append({"page": page_no, "source": "reaction",
                             "smiles": rxn.as_reaction_smiles(), "valid": ""})

            # Models miss structures silently, so every page gets a way to add them.
            if result.engine != "manual":
                st.caption("⚠️ Models can miss structures — compare with the "
                           "original page and add any that are missing.")
            added = st.text_area(
                "Add structures (one SMILES per line)", key=f"add_{page.index}",
                placeholder="e.g. CC(=O)Oc1ccccc1C(=O)O  (aspirin)",
            )
            n_added = 0
            for ai, smi in enumerate(added.split()):
                v = validate_smiles(smi)
                if not v.is_valid:
                    st.error(f"`{smi}` is not a valid SMILES — check it for typos.")
                    continue
                n_added += 1
                st.markdown(f"**➕ Added {n_added}** — `{v.canonical_smiles}`")
                # SMILES in the key: a changed line gets a fresh canvas.
                _review(v.canonical_smiles, f"canvas_add_{page.index}_{ai}_{smi}",
                        page_no, "added", rows)

    with export_slot:
        st.download_button(
            f"Export reviewed set ({len(rows)}) · CSV", _export_csv(rows),
            file_name="chemora_structures.csv", mime="text/csv",
            use_container_width=True, disabled=not rows,
        )


# ---------------------------------------------------------------------------
# Main flow
# ---------------------------------------------------------------------------
if run:
    if uploaded is None:
        st.error("Please upload a PDF first.")
    else:
        pdf_bytes = uploaded.read()
        progress = st.progress(0.0, text="Starting…")

        def _cb(frac, msg):
            progress.progress(min(max(frac, 0.0), 1.0), text=msg)

        with st.spinner("Processing… (first run loads models — this can take a while)"):
            result = process_pdf(pdf_bytes, progress_callback=_cb, engine=engine_choice)

        progress.empty()
        # Stash in session_state so edits/reruns don't reprocess the PDF.
        st.session_state["result"] = result

# Show the most recent result (persists across Streamlit reruns).
if "result" in st.session_state:
    _show_results(st.session_state["result"])
else:
    st.info(
        "👈 Upload a chemistry PDF and click **Run extraction** to begin. "
        "No GPU? Choose the **manual** engine to demo the canvas instantly."
    )
