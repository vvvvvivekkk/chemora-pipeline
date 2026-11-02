# Chemora Extract — phase-wise plan

Based on the 20-step brief from SciChem AI (Nov 2026). Each of the 20 steps is
mapped to a phase below, with what is already done in this repo and what is not.
Work the phases in order. Nothing in a later phase starts until the earlier
phase's "done when" is true, except the daily habits, which start now.

Glossary for beginners is at the bottom.

---

## Where we are (2026-11-02)

| Brief step | Status | Evidence in repo |
|---|---|---|
| 1–5 Understand first | **partly** | `docs/RESEARCH.md` covers OCSR; chemistry basics and SMILES practice not yet done by the developer himself |
| 6 One env per model, pinned | **no** | single `chemora` conda env; `requirements.txt` has loose pins (`>=`) |
| 7 No weights in git | **yes** | `.gitignore` ignores `*.ckpt *.pt *.pth models/ .cache/` |
| 8 Git properly (branch, small commits) | **partly** | small commits yes; everything goes straight to `main` |
| 9 Licence table | **no** | — |
| 10 Prove offline | **no** | weights download on first start; never tested with network off |
| 11 Benchmark 20+ images, verified keys | **partly** | 47 synthetic molecules (drawn from known SMILES) + 1 real paper (52 structures, 10 reactions) with OPSIN-derived key `bench/keys/bjoc-20-228.json`; no hand-drawn / phone-photo cases |
| 12 Benchmark script, InChIKey, time per image | **mostly** | `bench/benchmark.py`, `bench/evaluate.py` compare by InChIKey; time per image not reported |
| 13 "valid ≠ correct" understood | **yes** | README says so; UI still shows ✅ which can be misread as "correct" |
| 14 Failure log | **partly** | failures live in commit messages and `bench/scoreboard.md`; no dedicated log |
| 15 2–3 models + voting | **no** | DECIMER engine exists as an option but runs alone; no agreement logic |
| 16 MolDetect numbers (found / missed / false) | **partly** | evaluate.py reports false positives and missed on bjoc-20-228 only |
| 17 Conditions OCR, two engines, abbreviations | **no** | RxnScribe returns conditions text; it is never OCR'd or checked |
| 18 Reaction roles and labels | **partly** | reactants/products split exists; compound labels ("1", "4h") not linked to structures |
| 19 Save chemist corrections | **no** | edits are exported to CSV only; crop + model answers + reviewer not saved |
| 20 Daily numeric report, ask before core changes | **starting** | `docs/DAILY_LOG.md` created today |

Current score (commit `4cc4f43`, cycle 6): structures exact **44/52 (84.6%)**,
reactions exact **7/10**, false positives 0 structures / 1 reaction, synthetic
35/47 exact.

---

## Phase 0 — Housekeeping and understanding (this week)

Brief steps 1–5, 8, 14, 20.

Goal: the repo follows the owner's rules and the developer can explain the
pipeline without a tool open.

Tasks
- [x] `CLAUDE.md` with the five daily rules and the working rules.
- [x] `docs/DAILY_LOG.md` and `docs/FAILURE_LOG.md` started.
- [ ] Finish cycle 7 (substitute-label probe for R′/R″ pictures) that was
      running when the machine went offline; score it; keep or revert; log it.
- [ ] From now on: every change on a branch (`cycle-8-xxx`, `phase1-envs`),
      merged to `main` only with a scoreboard line that is not worse.
- [ ] Developer, by hand (no Claude): draw the pipeline on paper
      (input → cleanup → find → read → read text → assemble reaction →
      validate → review → save) and the 5 rules.
- [ ] Developer: with RDKit, turn 5 SMILES into pictures and 5 pictures (from
      the benchmark cases) back into SMILES. Script goes in `notebooks/` or
      `docs/learning/`.
- [ ] Developer: take one Suzuki coupling scheme from `bench/cases/` and label
      every part (reactant, reactant, arrow, conditions, product) in a note.

Done when: the four learning items are written in the developer's own words in
`docs/learning/`, and cycle 7 has a scoreboard line and a log entry.

---

## Phase 1 — Set up safely

Brief steps 6, 7, 9, 10.

Goal: the pipeline is reproducible, licence-checked and provably offline.

Tasks
- [ ] **Pin everything.** Freeze the working `chemora` env to
      `envs/openchemie/requirements.lock.txt` (exact `==` versions). Keep the
      human-readable `requirements.txt` as the explanation, the lock file as the
      truth.
- [ ] **One env per model family.** OpenChemIE already bundles MolDetect +
      MolScribe + RxnScribe on one PyTorch 1.13 stack, so that stays one env.
      Every *additional* model (DECIMER on TensorFlow, MolNexTR, an OCR engine)
      gets its own env under `envs/<name>/` with its own lock file. Models in
      different envs talk through a small subprocess call that passes image
      paths in and JSON out (`pipeline/engines/<name>_runner.py`).
- [ ] **Weights offline.** A one-time `scripts/fetch_weights.py` downloads every
      checkpoint into `models/` (gitignored) and prints each file's size and
      SHA-256. The pipeline loads from `models/` only. Set
      `HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1` in `run_local.bat` and the Colab
      notebook. Use `curl -L -C - --retry 5` for stalled downloads and check the
      size afterwards.
- [ ] **Offline proof.** Disable the network, run `test_e2e.py` and
      `bench/evaluate.py`. Any attempt to download is a bug; fix it. Record the
      run (date, command, result) in `docs/OFFLINE_PROOF.md`.
- [ ] **Licence table** `docs/LICENSES.md`: tool, code licence, weights licence,
      commercial use yes/no/unsure, source link. Known so far — RDKit: BSD-3
      (yes). Streamlit: Apache-2.0 (yes). Ketcher: Apache-2.0 (yes).
      OpenChemIE / MolScribe / RxnScribe code: MIT (yes); **their weights:
      unsure, check each Hugging Face model card.** PyMuPDF (`fitz`):
      **AGPL-3.0 or paid commercial licence — a problem for a product; plan to
      swap to pdfium (`pypdfium2`, BSD/Apache) or poppler via `pdf2image`.**
      poppler/pdftotext: GPL (fine as a separate process, check). DECIMER code:
      MIT; weights: check. OPSIN: Artistic-2.0 (only used to build answer keys,
      not shipped). Anything "unsure" is blocked from the product until resolved.

Done when: a fresh machine can be set up from the lock files alone, the offline
run is recorded, and the licence table has no "unsure" rows for anything the app
imports.

---

## Phase 2 — Measure before improving

Brief steps 11–14, 16.

Goal: a benchmark the owner trusts, covering the hard cases, with every failure
written down.

Tasks
- [ ] **Extend the key set to 20+ real images of mixed difficulty** under
      `bench/keys/`: stereo wedges and hashes, abbreviations (Boc, Ph, OTf,
      dppf, CO₂Et), a hand-drawn structure, a phone photo of a printed page,
      a scanned lab-notebook page, a patent drawing. Each key SMILES from a
      verified source only: OPSIN on the paper's IUPAC name, PubChem CID, or a
      ChemDraw export. Record the source per entry (`"source": "opsin:<name>"`).
      Extend the existing builder pattern (`bench/keys/build_*.py`).
- [ ] Add the other three BJOC papers already in `bench/papers/`
      (bjoc.20.230, .233, .234) as keyed papers, starting with .233 (stacked
      CF₂ failure mode).
- [ ] **Time per image** and total time in the score line.
- [ ] **MolDetect numbers**: for every keyed paper report found / missed /
      falsely found boxes (precision and recall), separately from reading
      accuracy, so detection and reading failures are not confused.
- [ ] **Failure log automation**: `evaluate.py` writes every non-exact crop to
      `docs/FAILURE_LOG.md`-compatible rows (image path, model SMILES, correct
      SMILES, category left blank for the developer to fill). The developer
      writes the "probable reason" by hand.
- [ ] Rename the UI badge: ✅ → "parses" / ⚠️ → "does not parse". Never "correct".

Done when: ≥ 20 keyed images across ≥ 6 difficulty types, score line includes
timing and detection precision/recall, and every non-exact case has a reason in
the failure log.

---

## Phase 3 — Voting across 2–3 structure readers

Brief step 15. **Ask the owner before starting; this changes core logic.**

Goal: confidence comes from agreement between independent models, not from
RDKit parsing.

Tasks
- [ ] Add DECIMER (TensorFlow) and MolNexTR as readers in their own envs
      (Phase 1 pattern). Each takes the same MolDetect crop and returns SMILES.
- [ ] `pipeline/vote.py`: compare by InChIKey. All agree → 🟢 trust; two of
      three agree → 🟡 check; all differ → 🔴 redraw. Separate flag
      "agree on skeleton, differ on stereo" (first 14 characters of the InChIKey
      match, full key differs).
- [ ] Post-processing (`postprocess_smiles`) runs on every reader's output
      before voting, so a fix helps all of them.
- [ ] Score each reader alone *and* the vote on the benchmark. Record whether
      "all agree" cases are actually correct against the key — the owner's
      warning: all models can agree and still be wrong.
- [ ] UI and CSV carry the 🟢/🟡/🔴 level and each reader's answer.

Done when: the scoreboard has a line per reader and one for the vote, and the
🟢 bucket's accuracy against the key is reported.

---

## Phase 4 — Conditions text OCR

Brief step 17.

Goal: the text above the arrow (reagents, solvent, temperature, yield) is read
and checked.

Tasks
- [ ] Crop the conditions region RxnScribe already locates; run two offline OCR
      engines (Tesseract and PaddleOCR or EasyOCR, each in its own env).
- [ ] Build a small keyed set of conditions crops (≥ 30) with the correct text
      typed from the paper's own text layer (PyMuPDF/pdftotext), which is
      ground truth, not a model.
- [ ] Compare the two engines; known trouble spots: subscripts (H₂O), °C, "Pd",
      Greek letters.
- [ ] `pipeline/chem_text.py`: abbreviation dictionary (Boc, Ph, OTf, dppf,
      DMF, THF, rt, equiv, …) and misread rules (e.g. `0C` → `°C`,
      `H20` → `H₂O`). Dictionary lives in a JSON file, not in code.
- [ ] Export conditions as structured fields: reagents, solvent, temperature,
      time, yield.

Done when: conditions exact-match rate per engine is on the scoreboard and the
better engine is wired in.

---

## Phase 5 — Reaction roles and compound labels

Brief step 18.

Goal: every structure knows its label ("1", "4h") and every reaction knows
which structure is which.

Tasks
- [ ] Detect compound labels near each box (OCR from Phase 4 on the strip
      below/beside the box; the PDF text layer when available).
- [ ] Link label ↔ structure ↔ reaction role (reactant / reagent / product).
- [ ] Add label accuracy to the bjoc-20-228 key (labels are already in the
      paper) and to the score line.
- [ ] CSV gains `label` and `role` columns.

Done when: labels are scored and ≥ 90% correct on the keyed papers.

---

## Phase 6 — Capture chemist corrections

Brief step 19.

Goal: every review action becomes future training data.

Tasks
- [ ] On "Apply" in the canvas, save to `data/corrections/<date>/`: the crop
      image, each reader's SMILES, the vote level, the corrected SMILES, the
      timestamp and the reviewer name. One JSON per correction plus the PNG.
- [ ] Same for "added" (model missed it) and "deleted" (false positive).
- [ ] A script that turns corrections into new benchmark key entries after a
      second reviewer confirms them.

Done when: a review session produces a corrections folder that
`bench/evaluate.py` can score against.

---

## Daily habits (start now)

- End of each session: a `docs/DAILY_LOG.md` entry in numbers — what was tried,
  score line before/after, new failures, blockers, what was learned (in the
  developer's own words).
- Every reverted attempt → `docs/FAILURE_LOG.md`.
- Before touching voting, validation, scoring or keys → ask the owner.
- Start every Claude session by reading `CLAUDE.md` (it holds the five rules).

---

## Glossary (plain words)

- **OCSR** — optical chemical structure recognition: picture of a molecule → text.
- **SMILES** — a text spelling of a molecule, e.g. `Fc1ccccc1Br` (a benzene ring
  with one F and one Br). Many spellings can mean the same molecule.
- **InChIKey** — a fixed-length fingerprint of a molecule. Same molecule → same
  key, whatever the SMILES spelling. First 14 characters = atoms and bonds
  only; the rest adds stereo and charge. We compare molecules with it.
- **Canonical SMILES** — RDKit's one standard spelling for a molecule.
- **Valid (RDKit)** — the SMILES can be parsed. Says nothing about whether it is
  the molecule in the picture.
- **Bounding box (bbox)** — the rectangle MolDetect draws around one structure.
- **Post-processing** — code that fixes known systematic model errors after
  reading, before validation (`postprocess_smiles` in `molecule_engine.py`).
- **Cache** — saved raw model outputs so post-processing changes can be
  re-scored in seconds without re-running the models (~13 min).
- **Valence** — how many bonds an atom can have (carbon 4). Over-valent =
  impossible, RDKit rejects it.
- **R-group** — a placeholder in a drawing ("some group goes here"); `*` in SMILES.
- **Reactant / product / reagent / solvent / yield** — what goes in, what comes
  out, what helps the reaction, what it is dissolved in, how much product was
  obtained (%).
- **Stereochemistry (wedge/hash bonds)** — bold wedge = bond coming toward you,
  hashed = going away; it distinguishes mirror-image molecules.
- **Precision / recall (detection)** — of the boxes found, how many were real
  (precision); of the real structures, how many were found (recall).
