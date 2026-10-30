"""Score the pipeline on a real paper against a hand-built answer key.

    python bench/evaluate.py bjoc-20-228

Reads the cached model output (bench/cache/<paper>/, filled once by the model),
runs the pipeline's post-processing on it, matches every crop to the answer key
(bench/keys/<paper>.json) by position, and writes:
  * bench/report.html  — original | model's drawing | correct drawing per crop,
                         reactions model vs correct, filter buttons;
  * one score line, printed and appended to bench/scoreboard.md with the commit.
The synthetic benchmark (also cached) is scored in the same run, so a change that
helps the paper but hurts the synthetic set shows up on the same line.
"""
from __future__ import annotations

import base64
import datetime
import html
import io
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

from benchmark import _key, _sim, model_outputs, run as run_synthetic  # noqa: E402

COLOR = {"exact": "#1a7f37", "near": "#b58100", "wrong": "#cf222e", "missed": "#cf222e",
         "fp": "#cf222e", "junk-dropped": "#6e7781", "invalid": "#cf222e"}


# ---------------------------------------------------------------------------
# grading
# ---------------------------------------------------------------------------
def grade(got: str | None, want: str) -> str:
    """exact = same molecule; near = same skeleton or Tanimoto >= 0.7; else wrong."""
    if not got:
        return "missed"
    kg, kw = _key(got), _key(want)
    if kg is None:
        return "invalid"
    if kg == kw:
        return "exact"
    if (len(kg) == 27 and len(kw) == 27 and kg[:14] == kw[:14]) or _sim(got, want) >= 0.7:
        return "near"
    return "wrong"


def _contains(box, x, y):
    return box[0] <= x <= box[2] and box[1] <= y <= box[3]


def _iou(a, b):
    ix = max(0.0, min(a[2], b[2]) - max(a[0], b[0]))
    iy = max(0.0, min(a[3], b[3]) - max(a[1], b[1]))
    inter = ix * iy
    union = (a[2] - a[0]) * (a[3] - a[1]) + (b[2] - b[0]) * (b[3] - b[1]) - inter
    return inter / union if union else 0.0


def match_structures(key, raw_mols, label_images=None):
    """Assign each crop to the key entry whose box holds the crop's centre (best IoU)."""
    from pipeline.molecule_engine import postprocess_smiles
    crops = []
    for p, items in enumerate(raw_mols):
        for i, m in enumerate(items):
            b = m["bbox"]
            cx, cy = (b[0] + b[2]) / 2, (b[1] + b[3]) / 2
            cands = [k for k in key["structures"] if k["page"] == p + 1 and _contains(k["box"], cx, cy)]
            best = max(cands, key=lambda k: _iou(k["box"], b), default=None)
            crops.append({"page": p + 1, "i": i, "bbox": b, "raw": m.get("smiles"),
                          "final": postprocess_smiles(m, label_images[p] if label_images else ()),
                          "crop": m.get("crop"),
                          "confidence": m.get("confidence"), "key": best})
    # one crop per key entry: the best-overlapping; other crops on it are duplicates
    taken = {}
    for c in sorted(crops, key=lambda c: -_iou(c["key"]["box"], c["bbox"]) if c["key"] else 0):
        k = c["key"]
        if k is None:
            c["status"] = "fp" if c["final"] else "junk-dropped"
            c["why"] = "matches nothing in the paper"
        elif k["junk"]:
            c["status"] = "fp" if c["final"] else "junk-dropped"
            c["why"] = f"junk region: {k['label']}"
        elif k["id"] in taken:
            c["status"] = "fp" if c["final"] else "junk-dropped"
            c["why"] = f"duplicate crop of {k['label']}"
        else:
            taken[k["id"]] = c
            c["status"] = grade(c["final"], k["smiles"])
            c["why"] = "dropped by cleanup" if c["raw"] and not c["final"] else ""
    rows = []
    for k in key["structures"]:
        if k["junk"]:
            continue
        c = taken.get(k["id"])
        rows.append({"key": k, "crop": c, "status": c["status"] if c else "missed"})
    fps = [c for c in crops if c["status"] in ("fp", "junk-dropped")]
    return rows, fps


def match_reactions(key, raw_rxns, label_images=None, page_sizes=None):
    from pipeline.reaction_engine import reactions_from_raw
    model = reactions_from_raw(raw_rxns, label_images, page_sizes)
    keyset = lambda smis: sorted(_key(s) or s for s in smis)
    used, rows = set(), []
    for k in key["reactions"]:
        want_r, want_p = keyset(k["reactants"]), keyset(k["products"])
        best, status = None, "missed"
        for j, r in enumerate(model):
            if j in used or r.page != k["page"] - 1:
                continue
            got_r, got_p = keyset(r.reactants), keyset(r.products)
            s = "exact" if (got_r, got_p) == (want_r, want_p) else \
                "near" if got_p == want_p or got_r == want_r else "wrong"
            rank = {"exact": 3, "near": 2, "wrong": 1}[s]
            if best is None or rank > best[0]:
                best = (rank, j, s)
        if best:
            used.add(best[1])
            status = best[2]
        rows.append({"key": k, "model": model[best[1]] if best else None, "status": status})
    fps = [r for j, r in enumerate(model) if j not in used]
    return rows, fps


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def _png(img, size=(260, 220)):
    if img is None:
        return ""
    img = img.copy()
    img.thumbnail(size)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def _mol_png(smiles):
    from rdkit import Chem
    from rdkit.Chem import Draw
    m = Chem.MolFromSmiles(smiles or "")
    return _png(Draw.MolToImage(m, size=(320, 260))) if m else ""


def _rxn_png(reactants, products):
    from rdkit.Chem import AllChem, Draw
    try:
        r = AllChem.ReactionFromSmarts(f"{'.'.join(reactants)}>>{'.'.join(products)}", useSmiles=True)
        return _png(Draw.ReactionToImage(r, subImgSize=(220, 200)), size=(900, 220))
    except Exception:
        return ""


def _cell(src, text=""):
    img = f'<img src="{src}">' if src else '<div class="none">—</div>'
    return f'<td>{img}<code>{html.escape(text or "")}</code></td>'


def write_report(paper, score, rows, fps, rrows, rfps, pages):
    def region(page, box):
        im = pages[page - 1].image
        W, H = im.size
        return im.crop((int(box[0] * W), int(box[1] * H), int(box[2] * W), int(box[3] * H)))

    out = []
    for r in rows:
        k, c, st = r["key"], r["crop"], r["status"]
        orig = _png(c["crop"]) if c and c["crop"] is not None else _png(region(k["page"], k["box"]))
        model_smi = (c["final"] or f"(dropped) {c['raw']}") if c else "(not detected)"
        note = c["why"] if c and c["why"] else ""
        out.append(f'<tr class="row {st}" data-kind="{"bad" if st not in ("exact",) else "ok"}">'
                   f'<td class="tag" style="background:{COLOR[st]}">{st}</td>'
                   f'<td>p{k["page"]} <b>{html.escape(k["label"])}</b><br><small>{html.escape(note)}</small></td>'
                   + _cell(orig) + _cell(_mol_png(c["final"]) if c else "", model_smi)
                   + _cell(_mol_png(k["smiles"]), k["smiles"]) + "</tr>")
    for c in fps:
        st = c["status"]
        out.append(f'<tr class="row {st}" data-kind="{"fp" if st == "fp" else "ok"}">'
                   f'<td class="tag" style="background:{COLOR[st]}">{"false positive" if st == "fp" else "junk dropped"}</td>'
                   f'<td>p{c["page"]}<br><small>{html.escape(c["why"])}</small></td>'
                   + _cell(_png(c["crop"])) + _cell(_mol_png(c["final"]), c["final"] or f"(dropped) {c['raw']}")
                   + _cell("", "nothing should be extracted here") + "</tr>")
    rx = []
    for r in rrows:
        k, m, st = r["key"], r["model"], r["status"]
        rx.append(f'<tr class="row {st}" data-kind="{"bad" if st != "exact" else "ok"}">'
                  f'<td class="tag" style="background:{COLOR[st]}">{st}</td>'
                  f'<td>p{k["page"]} <b>{k["id"]}</b></td>'
                  + _cell(_rxn_png(m.reactants, m.products) if m else "",
                          m.as_reaction_smiles() if m else "(not detected)")
                  + _cell(_rxn_png(k["reactants"], k["products"]),
                          f"{'.'.join(k['reactants'])}>>{'.'.join(k['products'])}") + "</tr>")
    for m in rfps:
        rx.append(f'<tr class="row fp" data-kind="fp"><td class="tag" style="background:{COLOR["fp"]}">false positive</td>'
                  f'<td>p{(m.page or 0) + 1}</td>' + _cell(_rxn_png(m.reactants, m.products), m.as_reaction_smiles())
                  + _cell("", "no such reaction in the paper") + "</tr>")
    (ROOT / "report.html").write_text(f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Extraction Report</title>
<style>
body{{font:14px system-ui,sans-serif;margin:16px;background:#fff;color:#1f2328}}
table{{border-collapse:collapse;width:100%}} td{{border-bottom:1px solid #d0d7de;padding:6px;vertical-align:top}}
td img{{max-width:260px;display:block}} code{{font-size:11px;word-break:break-all;display:block;max-width:300px}}
.tag{{color:#fff;font-weight:600;text-transform:uppercase;font-size:11px;width:70px}}
.score{{font:600 15px ui-monospace,monospace;background:#f6f8fa;padding:10px;border-radius:6px}}
button{{margin:8px 6px 8px 0;padding:6px 12px}} .none{{color:#8c959f}} th{{text-align:left}}
.hide{{display:none}}
</style></head><body>
<h2>{paper}: extraction vs answer key</h2>
<div class="score">{html.escape(score)}</div>
<div><button onclick="f('all')">All</button><button onclick="f('bad')">Wrong only</button>
<button onclick="f('fp')">False positives only</button></div>
<h3>Structures</h3>
<table><tr><th></th><th>Where</th><th>Original crop</th><th>Model's drawing (canvas)</th><th>Correct (answer key)</th></tr>
{''.join(out)}</table>
<h3>Reactions</h3>
<table><tr><th></th><th>Where</th><th>Model</th><th>Correct (answer key)</th></tr>{''.join(rx)}</table>
<script>
function f(k){{document.querySelectorAll('tr.row').forEach(r=>{{
 const d=r.dataset.kind; r.classList.toggle('hide', k==='all'?false: k==='bad'?!(d==='bad'||d==='fp'):d!=='fp');}});}}
</script></body></html>""", encoding="utf-8")


# ---------------------------------------------------------------------------
def _commit():
    try:
        h = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT.parent, text=True).strip()
        dirty = subprocess.check_output(["git", "status", "--porcelain", "--", "pipeline"], cwd=ROOT.parent, text=True).strip()
        return h + ("+dirty" if dirty else "")
    except Exception:
        return "no-git"


def main(paper):
    key = json.loads((ROOT / "keys" / f"{paper}.json").read_text())
    pages, raw_mols, raw_rxns = model_outputs(paper, (ROOT / "papers" / f"{paper}.pdf").read_bytes(),
                                              cache_only=True)
    rows, fps = match_structures(key, raw_mols, [pg.label_images for pg in pages])
    rrows, rfps = match_reactions(key, raw_rxns, [pg.label_images for pg in pages],
                                  [pg.image.size for pg in pages])
    syn = run_synthetic(cache_only=True)

    n = len(rows)
    count = lambda st: sum(r["status"] == st for r in rows)
    ex, rex = count("exact"), sum(r["status"] == "exact" for r in rrows)
    nfp = sum(c["status"] == "fp" for c in fps)
    score = (f"structures exact {ex}/{n} ({100 * ex / n:.1f}%), near {count('near')}, "
             f"wrong {count('wrong') + count('invalid')}, missed {count('missed')} · "
             f"reactions exact {rex}/{len(rrows)} ({100 * rex / len(rrows):.1f}%) · "
             f"false positives {nfp} structures, {len(rfps)} reactions · "
             f"synthetic exact {syn['exact']}/{syn['n']}, extra {syn['extra']}")
    write_report(paper, score, rows, fps, rrows, rfps, pages)
    line = f"| {datetime.date.today()} | `{_commit()}` | {paper} | {score} |"
    board = ROOT / "scoreboard.md"
    if not board.exists():
        board.write_text("# Scoreboard\n\n| date | commit | paper | score |\n|---|---|---|---|\n")
    with board.open("a", encoding="utf-8") as f:
        f.write(line + "\n")
    print("\n" + score)
    json.dump({r["key"]["id"]: r["status"] for r in rows} | {r["key"]["id"]: r["status"] for r in rrows},
              (ROOT / f"status_{paper}.json").open("w"), indent=0)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "bjoc-20-228")
