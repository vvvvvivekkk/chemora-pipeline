"""Ground-truth benchmark for the extraction pipeline.

We draw test pages from SMILES we already know, run them through the real
pipeline, and score the output against the answer key. Each case targets one
failure mode (tiny rings, stereo, charges, low resolution, reactions, ...).

    python bench/benchmark.py make     # build bench/cases/*.pdf + answers.json
    python bench/benchmark.py run      # run pipeline, write bench/results.json
    python bench/benchmark.py run small_rings stereo   # just some cases
    python bench/benchmark.py check    # self-check of the scoring logic
    python bench/benchmark.py papers   # real PDFs in bench/papers -> bench/papers_out

Run inside the `chemora` env from the project root.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
CASES_DIR = ROOT / "cases"
FONT = "C:/Windows/Fonts/arial.ttf"
PAGE = (2550, 3300)  # US letter at 300 DPI

LOREM = ("The title compounds were prepared as described in the general procedure "
         "and characterised by NMR and HRMS. Yields refer to isolated material. ")

# name -> list of pages. A page has molecules [(label, smiles)], optional
# reaction "A.B>>C", and render options: size (px), bw (black & white), text.
CASES = {
    "single_ibuprofen": [{"mols": [("Ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O")]}],
    "drug_grid": [{"mols": [
        ("Nicotine", "CN1CCC[C@H]1c1cccnc1"),
        ("Dopamine", "NCCc1ccc(O)c(O)c1"),
        ("Salicylic acid", "OC(=O)c1ccccc1O"),
        ("Fluoxetine", "CNCCC(Oc1ccc(cc1)C(F)(F)F)c1ccccc1"),
    ]}],
    "small_rings": [{"mols": [
        ("Benzene", "c1ccccc1"), ("Cyclohexane", "C1CCCCC1"),
        ("Pyridine", "c1ccncc1"), ("Naphthalene", "c1ccc2ccccc2c1"),
    ]}],
    "stereo": [{"mols": [
        ("L-Alanine", "N[C@@H](C)C(=O)O"),
        ("Menthol", "CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O"),
        ("Cholesterol", "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C"),
    ]}],
    "heteroatoms_charges": [{"mols": [
        ("Chloroquine", "CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12"),
        ("Betaine", "C[N+](C)(C)CC(=O)[O-]"),
        ("Sulfamethoxazole", "Cc1cc(no1)NS(=O)(=O)c1ccc(N)cc1"),
        ("Triphenylphosphine", "c1ccc(cc1)P(c1ccccc1)c1ccccc1"),
    ]}],
    "low_res_bw": [{"size": 220, "bw": True, "mols": [
        ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("Caffeine", "Cn1c(=O)c2c(ncn2C)n(C)c1=O"),
        ("Paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ]}],
    "paper_layout": [{"text": True, "mols": [
        ("1", "O=C1CCCCC1"), ("2", "OC1CCCCC1"),
    ]}],
    "rxn_esterification": [{"rxn": "OC(=O)c1ccccc1.CO>>COC(=O)c1ccccc1"}],
    "rxn_suzuki": [{"rxn": "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(cc1)-c1ccccc1"}],
    "multipage": [
        {"mols": [("Indole", "c1ccc2[nH]ccc2c1")]},
        {"mols": [("Quinoline", "c1ccc2ncccc2c1"), ("Phenol", "Oc1ccccc1")]},
    ],
    "text_only": [{"text": True, "mols": []}],
    "blank_then_mols": [  # an empty page must not wipe out the other pages
        {"text": True, "mols": []},
        {"mols": [("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"), ("Caffeine", "Cn1c(=O)c2c(ncn2C)n(C)c1=O")]},
    ],
    # --- drawing conventions: each must survive extraction AND the canvas ---
    "hashed_wedge_bonds": [{"mols": [  # dashed/hashed + solid wedges
        ("Glucose", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"),
        ("Ephedrine", "CN[C@@H](C)[C@H](O)c1ccccc1"),
    ]}],
    "explicit_h": [{"draw": "explicit_h", "mols": [
        ("Ethanol", "CCO"), ("Acetic acid", "CC(=O)O"),
    ]}],
    "abbreviations": [{"draw": "abbrev", "mols": [  # drawn as Ph / OMe / CO2Et
        ("Anisole-Ph", "COc1ccc(cc1)-c1ccccc1"),
        ("Ethyl ester", "CCOC(=O)c1ccc(OC)cc1"),
    ]}],
    "r_groups": [{"draw": "rgroup", "mols": [  # Markush: R1/R2 placeholders
        ("R-phenol", "*c1ccc(O)cc1"), ("R-amide", "*NC(=O)C*"),
    ]}],
    "dative_bond": [{"mols": [("BF3-NMe3", "C[N+](C)(C)[B-](F)(F)F")]}],
    # --- complex structures ---
    "complex_natural": [{"mols": [
        ("Morphine", "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5"),
        ("Strychnine", "O=C1C[C@H]2OCC=C3CN4CC[C@@]56[C@@H]4C[C@@H]3[C@@H]2[C@@H]6N1c1ccccc15"),
        ("Quinine", "COc1ccc2nccc([C@@H](O)[C@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1"),
        ("Camphor", "CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2"),
    ]}],
    "complex_rings": [{"mols": [
        ("Adamantane", "C1C2CC3CC1CC(C2)C3"),
        ("Spiro[4.5]decane", "C1CCC2(CC1)CCCC2"),
        ("Porphine", "c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3"),
        ("18-Crown-6", "C1COCCOCCOCCOCCOCCO1"),
    ]}],
    "complex_bio": [{"mols": [
        ("Gly-Ala-Phe", "NCC(=O)N[C@@H](C)C(=O)N[C@@H](Cc1ccccc1)C(=O)O"),
        ("Sucrose", "OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O"),
    ]}],
    "complex_paclitaxel": [{"size": 1300, "mols": [("Paclitaxel",
        "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C")]}],
    "metal_complex": [{"mols": [("Cisplatin", "[NH3]->[Pt](Cl)(Cl)<-[NH3]")]}],
}


# ---------------------------------------------------------------------------
# make: draw the pages
# ---------------------------------------------------------------------------
def _mol_image(smiles: str, size: int, bw: bool, draw: str | None = None):
    from rdkit import Chem
    from rdkit.Chem import Draw, rdAbbreviations
    opts = Draw.MolDrawOptions()
    if bw:
        opts.useBWAtomPalette()
    mol = Chem.MolFromSmiles(smiles)
    if draw == "explicit_h":
        mol = Chem.AddHs(mol)
    elif draw == "abbrev":
        mol = rdAbbreviations.CondenseMolAbbreviations(
            mol, rdAbbreviations.GetDefaultAbbreviations(), maxCoverage=0.6)
    elif draw == "rgroup":
        n = 0
        for a in mol.GetAtoms():
            if a.GetAtomicNum() == 0:
                n += 1
                a.SetProp("atomLabel", f"R{n}")
    return Draw.MolToImage(mol, size=(size, size), options=opts)


def _rxn_image(rxn: str):
    from rdkit.Chem import AllChem, Draw
    r = AllChem.ReactionFromSmarts(rxn, useSmiles=True)
    return Draw.ReactionToImage(r, subImgSize=(450, 450))


def _page(spec: dict):
    from PIL import Image, ImageDraw, ImageFont
    page = Image.new("RGB", PAGE, "white")
    d = ImageDraw.Draw(page)
    font, small = ImageFont.truetype(FONT, 42), ImageFont.truetype(FONT, 34)
    y = 200
    d.text((200, y), "Figure 1. Test structures.", fill="black", font=font)
    y += 120
    if spec.get("text"):
        for _ in range(8):
            d.text((200, y), LOREM[:95], fill="black", font=small)
            y += 55
        y += 60
    size = spec.get("size", 700)
    for i, (label, smi) in enumerate(spec.get("mols", [])):
        x = 200 + (i % 2) * 1100
        row_y = y + (i // 2) * (size + 200)
        page.paste(_mol_image(smi, size, spec.get("bw", False), spec.get("draw")), (x, row_y))
        d.text((x + 40, row_y + size + 20), label, fill="black", font=small)
    if spec.get("rxn"):
        img = _rxn_image(spec["rxn"])
        page.paste(img, (150, y + 200))
    if spec.get("text"):
        y = 2600
        for _ in range(6):
            d.text((200, y), LOREM[:95], fill="black", font=small)
            y += 55
    return page


def make():
    CASES_DIR.mkdir(exist_ok=True)
    for name, pages in CASES.items():
        imgs = [_page(p) for p in pages]
        imgs[0].save(CASES_DIR / f"{name}.pdf", "PDF", resolution=300,
                     save_all=True, append_images=imgs[1:])
    print(f"made {len(CASES)} cases in {CASES_DIR}")


# ---------------------------------------------------------------------------
# scoring
# ---------------------------------------------------------------------------
def _key(smiles: str) -> str | None:
    """InChIKey: the same molecule always gets the same key, however it's spelled."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles or "")
    if not m:
        return None
    # InChI can't encode R-groups (*): fall back to canonical SMILES for those.
    return Chem.MolToInchiKey(m) or Chem.MolToSmiles(m)


def _sim(a: str, b: str) -> float:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
    ma, mb = Chem.MolFromSmiles(a or ""), Chem.MolFromSmiles(b or "")
    if not (ma and mb):
        return 0.0
    fa, fb = (AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048) for m in (ma, mb))
    return DataStructs.TanimotoSimilarity(fa, fb)


def score_page(expected: list[str], predicted: list[str]) -> dict:
    """Match each expected molecule to the best unused prediction.

    exact     = same InChIKey (identical structure incl. stereo/charge)
    skeleton  = same connectivity, differs only in stereo/charge/isotope
    near      = Tanimoto >= 0.7 (recognisably the right molecule, with errors)
    missed    = nothing close
    extra     = predictions that matched nothing (hallucinations/duplicates)
    """
    unused = list(predicted)
    rows = []
    for exp in expected:
        ek = _key(exp)
        best, kind = None, "missed"
        for p in unused:
            if _key(p) == ek:
                best, kind = p, "exact"
                break
        if best is None:
            for p in unused:
                pk = _key(p)
                if pk and ek and pk[:14] == ek[:14]:
                    best, kind = p, "skeleton"
                    break
        if best is None and unused:
            cand = max(unused, key=lambda p: _sim(exp, p))
            if _sim(exp, cand) >= 0.7:
                best, kind = cand, "near"
        if best is not None:
            unused.remove(best)
        rows.append({"expected": exp, "got": best, "match": kind,
                     "sim": round(_sim(exp, best), 2) if best else 0.0})
    return {"rows": rows, "extra": unused}


def score_reaction(expected: str, predicted: list[tuple[list[str], list[str]]]) -> str:
    want_r, want_p = (sorted(_key(s) for s in side.split(".")) for side in expected.split(">>"))
    for reactants, products in predicted:
        got_r = sorted(_key(s) for s in reactants)
        got_p = sorted(_key(s) for s in products)
        if got_p == want_p and got_r == want_r:
            return "exact"
        if got_p == want_p:
            return "product_only"
    return "missed" if not predicted else "wrong"


def check():
    s = score_page(["c1ccccc1", "CCO", "N[C@@H](C)C(=O)O"],
                   ["C1=CC=CC=C1", "N[C@H](C)C(=O)O", "CCCCCC"])
    kinds = [r["match"] for r in s["rows"]]
    assert kinds == ["exact", "missed", "skeleton"], kinds
    assert s["extra"] == ["CCCCCC"], s["extra"]
    assert score_reaction("CO.OC=O>>COC=O", [(["OC=O", "CO"], ["COC=O"])]) == "exact"
    assert score_reaction("CO.OC=O>>COC=O", [(["CC"], ["COC=O"])]) == "product_only"
    assert score_reaction("CO>>CC", []) == "missed"
    assert score_page(["*c1ccccc1"], ["c1ccc(*)cc1"])["rows"][0]["match"] == "exact"
    from pipeline.molecule_engine import clean_model_smiles as clean
    assert clean("*.*") is None                       # logo
    assert clean("[Co].[Co].[Co].[Co]") is None       # licence badge
    assert clean("[Fl]c1ccccc1") == "*c1ccccc1"       # R' misread as flerovium
    assert clean("[Pa]c1ccccc1") == "*c1ccccc1"       # ... or protactinium
    assert clean("Brc1ccccc1") == "Brc1ccccc1"        # real halogen untouched
    assert clean("CC(=O)O.[W]") == "CC(=O)O"          # stray label atom
    assert clean("CC(=O)[O-].[Na+]") == "CC(=O)[O-].[Na+]"  # real salt kept
    assert clean("C1CC(") == "C1CC("                  # invalid passes through to be flagged
    assert clean("CC.c1ccc(P(c2ccccc2)c2ccccc2)cc1") == "c1ccc(P(c2ccccc2)c2ccccc2)cc1"  # neighbour leak
    print("scoring self-check OK")


# ---------------------------------------------------------------------------
# run: pipeline vs answer key
# ---------------------------------------------------------------------------
def run(engine: str = "openchemie", only: list[str] | None = None):
    from pipeline.orchestrator import process_pdf
    out, t_all = [], time.time()
    for name, pages in CASES.items():
        if only and name not in only:
            continue
        t0 = time.time()
        res = process_pdf((CASES_DIR / f"{name}.pdf").read_bytes(), engine=engine)
        dt = time.time() - t0
        case = {"case": name, "seconds": round(dt, 1), "status": res.status, "pages": []}
        for i, spec in enumerate(pages):
            pr = res.pages[i] if i < len(res.pages) else None
            got = [m.canonical_smiles or m.raw_smiles for m in pr.molecules] if pr else []
            page = {"page": i + 1, **score_page([s for _, s in spec.get("mols", [])], got)}
            if spec.get("rxn"):
                rx = [(r.reactants, r.products) for r in pr.reactions] if pr else []
                page["reaction"] = {"expected": spec["rxn"],
                                    "got": [f"{'.'.join(a)}>>{'.'.join(b)}" for a, b in rx],
                                    "match": score_reaction(spec["rxn"], rx)}
            case["pages"].append(page)
        out.append(case)
        print(f"{name:22s} {dt:6.1f}s  " + "  ".join(
            f"p{p['page']}:" + ",".join(r["match"] for r in p["rows"])
            + (f" extra={len(p['extra'])}" if p["extra"] else "")
            + (f" rxn={p['reaction']['match']}" if "reaction" in p else "")
            for p in case["pages"]), flush=True)
    (ROOT / "results.json").write_text(json.dumps(out, indent=2))

    rows = [r for c in out for p in c["pages"] for r in p["rows"]]
    n = len(rows)
    count = lambda k: sum(r["match"] == k for r in rows)
    extra = sum(len(p["extra"]) for c in out for p in c["pages"])
    print(f"\n{n} molecules: exact {count('exact')}  skeleton {count('skeleton')}  "
          f"near {count('near')}  missed {count('missed')}  |  extra {extra}  "
          f"|  total {time.time() - t_all:.0f}s")


# ---------------------------------------------------------------------------
# papers: real PDFs have no answer key, so we "reverse-engineer" each result:
# redraw what the model read and put it next to the original crop. A human (or
# a second model) can then check every structure at a glance.
# ---------------------------------------------------------------------------
def papers(only: list[str] | None = None):
    from PIL import Image, ImageDraw, ImageFont
    from pipeline import pdf_utils
    from pipeline.molecule_engine import extract_molecules_from_images
    from pipeline.reaction_engine import extract_reactions_from_images
    from pipeline.validate import validate_smiles

    out_dir = ROOT / "papers_out"
    out_dir.mkdir(exist_ok=True)
    font = ImageFont.truetype(FONT, 16)
    summary = []
    for pdf in sorted((ROOT / "papers").glob("*.pdf")):
        if only and pdf.stem not in only:
            continue
        t0 = time.time()
        pages = pdf_utils.pdf_to_images(pdf.read_bytes())
        rows, tiles, rxns = [], [], []
        for pg in pages:  # one page at a time: this machine has little RAM
            for d in extract_molecules_from_images([pg.image]):
                v = validate_smiles(d.smiles)
                rows.append({"page": pg.index + 1, "smiles": v.canonical_smiles or d.smiles,
                             "valid": v.is_valid, "confidence": round(d.confidence, 3),
                             "has_rgroup": "*" in d.smiles})
                tile = Image.new("RGB", (620, 330), "white")
                if d.crop is not None:
                    c = d.crop.copy(); c.thumbnail((300, 300)); tile.paste(c, (5, 25))
                if v.image is not None:
                    r = v.image.copy(); r.thumbnail((300, 300)); tile.paste(r, (315, 25))
                ImageDraw.Draw(tile).text(
                    (5, 3), f"p{pg.index + 1}  conf {d.confidence:.2f}  "
                    + ("" if v.is_valid else "INVALID  ") + (d.smiles[:48]),
                    fill="black" if v.is_valid else "red", font=font)
                tiles.append(tile)
            for r in extract_reactions_from_images([pg.image]):
                rxns.append({"page": pg.index + 1, "reaction": r.as_reaction_smiles(),
                             "conditions": r.conditions})
        for k in range(0, len(tiles), 12):  # contact sheets, 12 structures each
            sheet = Image.new("RGB", (1240, 330 * 6), "#ddd")
            for i, t in enumerate(tiles[k:k + 12]):
                sheet.paste(t, ((i % 2) * 620, (i // 2) * 330))
            sheet.save(out_dir / f"{pdf.stem}_{k // 12 + 1}.png")
        valid = sum(r["valid"] for r in rows)
        summary.append({"paper": pdf.name, "pages": len(pages), "molecules": len(rows),
                        "valid": valid, "low_conf": sum(r["confidence"] < 0.5 for r in rows),
                        "rgroup": sum(r["has_rgroup"] for r in rows),
                        "reactions": len(rxns), "seconds": round(time.time() - t0)})
        (out_dir / f"{pdf.stem}.json").write_text(json.dumps({"molecules": rows, "reactions": rxns}, indent=2))
        print(summary[-1], flush=True)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "run":
        run(only=sys.argv[2:] or None)  # e.g.  run small_rings stereo
    elif cmd == "papers":
        papers(only=sys.argv[2:] or None)  # e.g.  papers bjoc-20-228
    else:
        {"make": make, "check": check}[cmd]()
