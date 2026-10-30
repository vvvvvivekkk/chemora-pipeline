---
title: Chemora Extract
emoji: 🧪
colorFrom: indigo
colorTo: green
sdk: streamlit
sdk_version: 1.39.0
app_file: app.py
pinned: false
---

# Chemora Extract — PDF → Chemical Structures & Reactions

Reads a chemistry PDF (paper, patent, lab report), finds the chemical structures
and reaction schemes inside it, converts them to machine-readable SMILES, and
puts each one on an **editable canvas** where a chemist confirms or fixes it,
then exports the reviewed set as CSV.

> **Design principle — human-in-the-loop.** No optical structure recognition
> (OCSR) model is 100% accurate, not even commercial ones (see the measured
> numbers below). The AI does the heavy lifting; the chemist reviews on the
> canvas. The canvas is the correction step, not just the output.

---

## Quick start (Windows)

Needs [Anaconda/Miniconda](https://www.anaconda.com/download) and **≥ 8 GB RAM**
(the three models take ~1–2 GB; close other heavy apps).

```bat
git clone https://github.com/vvvvvivekkk/chemora-pipeline.git
cd chemora-pipeline
run_local.bat
```

The first run creates a `chemora` conda env (Python 3.9, PyTorch 1.13 CPU,
OpenChemIE, RDKit — takes a while), then opens the app at
<http://localhost:8501>. Later runs just activate the env and launch.

**Manual setup** (any OS — this is exactly what `run_local.bat` does):

```bash
conda create -y -n chemora python=3.9
conda activate chemora
conda install -y -c conda-forge poppler pdftotext rdkit
pip install torch==1.13.1 --index-url https://download.pytorch.org/whl/cpu
pip install "OpenChemIE @ git+https://github.com/CrystalEye42/OpenChemIE"
pip install pymupdf streamlit streamlit-ketcher
streamlit run app.py
```

Install order matters: torch before OpenChemIE (it needs torch < 2), and
OpenChemIE before the app deps (so its Pillow/OpenCV pins win).

**No model / instant demo:** pick **manual** in the engine dropdown (or set
`CHEMORA_ENGINE=manual`). No model loads; you type SMILES and edit them on the
canvas. Needs only the core packages.

### Settings (environment variables)

| Variable | Default | Meaning |
|---|---|---|
| `CHEMORA_ENGINE` | `openchemie` | `openchemie` (molecules + reactions), `decimer` (molecules only, needs TensorFlow), `manual` |
| `CHEMORA_MAX_PAGES` | `10` | Only the first N pages are processed; `0` = all |
| `CHEMORA_PDF_DPI` | `300` | Page render resolution |
| `CHEMORA_DEVICE` | `auto` | `cuda`, `cpu` or `auto` |
| `CHEMORA_LOW_CONF` | `0.5` | Below this model confidence a result is flagged |
| `CHEMORA_CROP_PAD` | `0` | Grow each detected box by this fraction before reading it (0.08 tested: no gain on real papers, pulled in captions) |

---

## Using the app

1. **Upload** a PDF in the sidebar and click **Run extraction** (first run
   loads the models: ~1 min on CPU, then ~15–30 s per page).
2. Each page shows the **original on the left** and the extracted structures
   on the right — ✅ valid, ⚠️ RDKit could not parse it.
3. **Fix a structure** on its canvas, then click **Apply**. The app shows
   `✏️ Edited → <new SMILES>`. **Reset** restores the model's version.
4. **Missed one?** Every page has an *Add structures* box (one SMILES per line).
   Typos are rejected with a message instead of a blank canvas.
5. **Export reviewed set · CSV** (sidebar) downloads what you ended up with:
   `page, source (extracted / extracted+edited / added / reaction), smiles, valid`.

---

## How it works

```
PDF
 │  pdf_utils.py         PyMuPDF renders each page to a 300-DPI image
 ▼
page images
 │  molecule_engine.py   MolDetect  — finds a box around every structure
 │                       MolScribe  — reads each box: image → atom/bond graph
 │                                    → SMILES + molfile + confidence
 │                       cleanup    — fixes systematic misreads (below)
 │  reaction_engine.py   RxnScribe  — finds reaction schemes:
 │                                    reactants >> products + conditions text
 ▼
raw SMILES
 │  validate.py          RDKit — is it real chemistry? canonical form + picture
 ▼
results
 │  orchestrator.py      the one function the UI calls: process_pdf()
 │  app.py               Streamlit: page | Ketcher canvas per structure,
 ▼                       edits + additions tracked, CSV export
reviewed set (CSV)
```

* **MolDetect, MolScribe, RxnScribe** are the three models inside
  [OpenChemIE](https://github.com/CrystalEye42/OpenChemIE) (MIT, Coley &
  Barzilay groups). Everything runs locally — no data leaves the machine.
* MolScribe predicts an atom/bond **graph** (not just a string), so stereo
  wedges/hashes, charges and abbreviations (Ph, OMe, CO₂Et) are resolved into
  real atoms.
* The engine choice is passed per run (`process_pdf(..., engine=...)`), never
  written to shared config, so users of a shared server don't affect each other.
* A page with no structures returns nothing (MolScribe itself crashes on an
  empty batch; we short-circuit before calling it).
* **Cleanup of model output** (`clean_model_smiles`), from errors seen on real
  papers: R-group labels in coloured circles/boxes (R′, R″) come back as
  radioactive elements (Rn, Pa, Rf, Fl) → turned back into R-groups `*`; a crop
  is one molecule, so stray neutral fragments (a `Pd` catalyst label, `CH4`, the
  edge of a neighbouring drawing) are dropped while charged ions stay; anything
  with no carbon (journal logo, licence badge read as `Co.Co.Co.Co`) is dropped.

---

## Testing

### 1. End-to-end smoke test

```bash
python test_e2e.py          # sample PDF → pipeline → at least one valid molecule
```

### 2. Ground-truth benchmark (`bench/benchmark.py`)

Real papers don't come with an answer key, so we **reverse-engineer one**: draw
test pages from SMILES we already know, run them through the real pipeline,
and score the output against the SMILES we started from.

```bash
python bench/benchmark.py check            # self-test of the scoring logic
python bench/benchmark.py make             # bench/cases/*.pdf from known SMILES
python bench/benchmark.py run              # pipeline vs answer key → bench/results.json
python bench/benchmark.py run small_rings stereo   # only some cases
```

Scoring compares **InChIKeys** (the same molecule always gets the same key,
however the SMILES is spelled):
`exact` = identical incl. stereo/charge · `skeleton` = right atoms and bonds,
wrong stereo/charge · `near` = Tanimoto ≥ 0.7 · `missed` · `extra` = output
that matches nothing in the key.

Each case targets one failure mode: single molecules, grids, small rings,
stereo (wedge **and hashed/dashed** bonds), heteroatoms and charges, low-res
black-and-white, text around figures, explicit H, abbreviations, R-groups,
dative bonds, metal complexes, reactions, multi-page and empty pages, and
complex structures (morphine, strychnine, quinine, camphor, adamantane, spiro
rings, porphine, 18-crown-6, a tripeptide, sucrose, paclitaxel).

### 3. Real papers, cross-checked (`benchmark.py papers`)

```bash
# put PDFs in bench/papers/, then:
python bench/benchmark.py papers           # → bench/papers_out/
```

For every structure found, it writes the **original crop next to an RDKit
redraw of what the model read** (contact sheets, 12 per image) plus a JSON of
SMILES, validity, confidence and reactions. Mismatches are visible at a glance —
that's the cross-check when no answer key exists.

---

## Measured accuracy

### Synthetic benchmark (47 molecules, 22 cases, CPU)

**35 / 47 exact (74 %)**, 1 skeleton-only, 11 missed. Reactions: Suzuki exact,
esterification product-only.

| Case | Exact | | Case | Exact |
|---|---|---|---|---|
| single drug, drug grid | 5/5 | | hashed/wedge bonds | 1/2 |
| heteroatoms + charges | 4/4 | | small single rings | 1/4 |
| low-res black & white | 3/3 | | simple stereo | 2/3 |
| text around figure | 2/2 | | complex natural products | 3/4 |
| multi-page / blank page | 5/5 | | cages, macrocycles, crown | 1/4 |
| explicit H, abbreviations, R-groups | 6/6 | | paclitaxel (11 stereocentres) | 0/1 |
| tripeptide, sucrose | 2/2 | | metal / dative complexes | 0/2 |

### Real papers

Four open-access 2024 papers from the *Beilstein J. Org. Chem.* (CC-BY):
benzo[f]quinazolinediones (bjoc.20.228), spiroindolenines (bjoc.20.230),
chiral fluorinated amines (bjoc.20.233), porphyrin coupling (bjoc.20.234).

| Paper | Pages | Structures found | RDKit-valid | Low confidence | Reactions |
|---|---|---|---|---|---|
| bjoc.20.228 | 12 | 52 | 48 | 5 | 11 |
| bjoc.20.230 | 10 | 78 | 68 | 7 | 29 |
| bjoc.20.233 | 8 | 51 | 43 | 6 | 19 |
| bjoc.20.234 (porphyrins) | 15 | 63 | 43 | 20 | 8 |

"Valid" only means RDKit can parse it — **not** that it is the right molecule.
Correctness was graded by eye from the crop-vs-redraw sheets:

* **bjoc.20.228, pages 2–4: 11 / 24 correct (46 %).** On the same ten
  structures before the cleanup step: 4/10 + 2 false positives (journal logo,
  licence badge) → after: 6/10, no false positives.
* **bjoc.20.233 (one sheet, before cleanup): ~4 / 12.** This paper draws CF₂ as
  a stacked `C` with `F₂` above/below; MolScribe drops the fluorines
  (CF₂CF₂Br → CH₂CH₂Br). It also marks new stereocentres with a printed `*`,
  which is read as an R-group.

**Bottom line:** ~74 % exact on clean drawings, ~45 % on real journal figures
with coloured highlighting and publisher drawing styles. 98 % automatic accuracy
is not reachable with this model; the review canvas is what gets the final set
right. Ways to go higher: fine-tune MolScribe on these publisher styles
(coloured R-group circles, stacked CF₂, NMe₂), add a second engine (DECIMER)
and flag disagreements for review, or compare against a commercial engine.

### Known limitations (found by the benchmark)

* **Lone small rings** (benzene, cyclohexane, pyridine) are often not detected
  at all — MolDetect skips them. Use *Add structures* for these.
* **Terminal methyl on a wedge/hash bond** can be misread (as an alkyne or an
  R-group). Check stereocentres on the canvas.
* **Cages, macrocycles, crown ethers, very large molecules** (porphine,
  18-crown-6, paclitaxel) fail or come back garbled.
* **Metal complexes and dative bonds** are not detected.
* **R-group numbers** (R1, R2) are lost: SMILES only has `*`, so the canvas
  shows a generic `A` atom.
* **Coloured R-group circles/boxes** (R′, R″) are often read as a bonded Pd or
  deuterium; **NMe₂** and labels at the top edge of a highlighted ring are
  often read as an R-group.
* **Stacked labels** (C with F₂ written above/below) lose the stacked part.
* For reference, OpenChemIE's own paper reports F1 69.5 % on reaction schemes
  with R-groups and 64.3 % agreement with Reaxys
  ([arXiv:2404.01462](https://arxiv.org/abs/2404.01462)).

---

## Deploy

* **Hugging Face Spaces** (free CPU tier works, GPU faster): create a
  *Streamlit* Space, push this repo; it builds from `requirements.txt` +
  `packages.txt` and runs `app.py` (the YAML header at the top of this file is
  the Space config). Model weights download on first start.
* **Google Colab** (free GPU): see `notebooks/chemora_colab.ipynb`.
* **Any Linux server**: the manual setup above, then
  `streamlit run app.py --server.port 8501 --server.headless true`.

---

## Similar free tools online

| Tool | What it does |
|---|---|
| [DECIMER.ai](https://decimer.ai) | Upload a PDF or image → segmented structures → SMILES (Jena) |
| [MolScribe demo](https://huggingface.co/spaces/yujieq/MolScribe) | Single structure image → SMILES + molfile (the model this app uses) |
| [SciLifeLab Image-to-SMILES](https://image-to-smiles.serve.scilifelab.se) | Image or webcam photo → SMILES (DECIMER-based) |
| [OSRA](https://cactus.nci.nih.gov/osra/) | Classic rule-based image → SMILES (NIH) |
| [Mathpix](https://www.mathpix.com/image-to-smiles) | Screenshot → SMILES (commercial, has a free tier) |
| [Ketcher](https://lifescience.opensource.epam.com/ketcher/) | The open-source structure editor used for our canvas |
| Source: [OpenChemIE](https://github.com/CrystalEye42/OpenChemIE), [MolScribe](https://github.com/thomas0809/MolScribe), [DECIMER](https://github.com/Kohulan/DECIMER-Image_Transformer) | Run them yourself |

---

## File map

```
app.py                  Streamlit UI (upload, run, review on canvas, export)
run_local.bat           One-click Windows setup + launch
requirements.txt        Python dependencies (see install-order note inside)
packages.txt            System packages for HF Spaces (poppler, etc.)
pipeline/
  config.py             Engine + tunable settings (env vars above)
  pdf_utils.py          PDF → page images
  molecule_engine.py    Page image → SMILES (OpenChemIE / DECIMER)
  reaction_engine.py    Reaction schemes → reactants >> products (OpenChemIE)
  validate.py           RDKit: validate, canonicalize, draw
  orchestrator.py       Ties the stages together; the one API the UI calls
bench/benchmark.py      Ground-truth benchmark + real-paper cross-check
test_e2e.py             End-to-end smoke test
sample_data/            Sample PDF + where to get test papers
docs/RESEARCH.md        Background research
```
