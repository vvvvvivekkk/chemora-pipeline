# Research: PDF → Structure extraction, existing products, and pharma tooling

This document backs the Chemora Extract prototype. It covers (1) the technical
problem, (2) the open-source tools we build on, (3) the commercial products that
already do pieces of this, (4) what the pharma industry actually uses day to day,
(5) UI/UX patterns worth borrowing, and (6) where Chemora's real opportunity is.

> **One term up front — OCSR.** "Optical Chemical Structure Recognition" is the
> chemistry version of OCR: turning a *picture* of a molecule into a
> machine-readable structure (a SMILES string). "Document-to-structure" is OCSR
> applied to whole documents (PDFs, patents).

---

## 1. The problem, honestly

Extracting structures and reactions from a PDF is a real, hard research field —
not a weekend API hookup. The pipeline everyone uses is the same shape:

```
PDF → find the figures → segment individual structures/reactions →
recognise each (image→SMILES) → validate chemistry → (optionally) place on a canvas
```

**Accuracy is the crux.** The best published end-to-end system (OpenChemIE)
reports ~64% agreement against Reaxys' professionally curated extractions. No tool
— open or commercial — hits "zero data loss" today. **This is why the editable
canvas matters: the chemist confirms/fixes, and the AI just removes 90% of the
manual work.** Any pitch that promises 100% automated accuracy is a pitch that
will fail in a demo with a real chemist.

---

## 2. Open-source tools (what this prototype builds on)

| Tool | Does | Framework | Notes |
|---|---|---|---|
| **OpenChemIE** (MIT) | Molecules **and** reactions from a PDF | PyTorch | The closest single tool to Chemora's vision. Bundles MolScribe + RxnScribe + text parsing. Reads a PDF directly. **Our default engine.** |
| **DECIMER / DECIMER.ai** (Steinbeck lab) | PDF→segment→molecule SMILES | TensorFlow | Mature, well-documented, free web platform at decimer.ai. Molecule-only. Our fallback engine. |
| **MolScribe** (MIT) | Single molecule image → structure graph → SMILES | PyTorch | Image-to-graph; output is inspectable/correctable — fits human-in-the-loop. |
| **RxnScribe** (MIT) | Reaction diagram → reactants/conditions/products | PyTorch | The reaction-parsing piece inside OpenChemIE. |
| **MolNexTR, MolGrapher, MolScribe** | Newer OCSR models | PyTorch | Active research frontier; accuracy climbing. |
| **OSRA, MolVec, Imago** | Classic rule-based OCSR | C++/Java | Older; generally beaten by the deep-learning models on messy real figures, but lightweight. |
| **RDKit** | Validate / canonicalize / draw SMILES | C++/Python | The standard open-source chemistry engine. Our validation layer. |
| **Ketcher** (EPAM, Apache-2.0) | Web chemical drawing canvas | TypeScript/React | Imports/exports SMILES, MOL, InChI, RXN; has a React component + a `getSmiles()/setMolecule()` API. **This is almost certainly what CanvasPrime is built on** (identical tool-rail + element-sidebar + ring-template layout). |

**Why this matters for the role:** Chemora doesn't need to invent OCSR. The
engineering job is to *assemble, serve, and productize* these — the exact
"familiar with these models, APIs, deployment, client handover" ask.

---

## 3. Commercial products (what Chemora competes with / positions against)

| Product | Vendor | What it does | Why it's not enough for the target market |
|---|---|---|---|
| **Reaxys** | Elsevier | Massive curated DB; extracts substances/reactions from literature & patents; AI retrosynthesis on 73M+ reactions; commercial availability + pricing of reagents | Expensive enterprise subscription; cloud-only; a *database* product, not "drop my own PDF → editable canvas." Out of reach for mid/small Indian pharma. |
| **SciFinder-n** | CAS (ACS) | The other giant literature/substance/reaction search; retrosynthesis | Same: premium subscription, cloud, not a self-serve extraction-to-canvas tool. |
| **ChemDraw** (Signals ChemDraw) | Revvity (ex-PerkinElmer) | The *industry-standard* drawing tool; "Structure from Name," some image→structure | Per-seat licensing; drawing-first, not a document-extraction engine; pricey at scale. |
| **MarvinSketch / "Document to Structure"** | ChemAxon | Drawing + a document-to-structure OCSR feature | Enterprise licensing; part of a larger JChem stack. |
| **CLiDE / CLiDE Pro** | Keymodule | Dedicated commercial OCSR (document→structure) | Niche, licensed; not a modern self-serve web product. |
| **IBM RXN for Chemistry** | IBM | Free-ish web tool: reaction prediction + retrosynthesis, clean UX | Not focused on PDF structure extraction; cloud; not on-prem. |

**Takeaway:** The capability exists at the top end (Reaxys/SciFinder) and in the
drawing tools (ChemDraw/Marvin), but it's **expensive, cloud-based, and
enterprise-shaped**. Nobody is serving the "affordable, self-serve, on-prem,
PDF→editable-canvas" slice for mid-tier labs — which is exactly the founder's
thesis.

---

## 4. What the pharma industry actually uses day to day

A working medicinal/process chemist's stack, roughly by category:

- **Drawing structures:** **ChemDraw** (dominant, almost universal in pharma),
  **MarvinSketch** (ChemAxon), **Ketcher** (open, increasingly embedded in web
  tools and ELNs).
- **Searching literature / reactions / substances:** **Reaxys** (Elsevier),
  **SciFinder-n** (CAS), **PubChem** (free), **Google Scholar / patents**.
- **Electronic Lab Notebooks (ELN) & data management:** **Benchling** (modern,
  popular with biotech/startups), **IDBS E-WorkBook**, **Revvity Signals**,
  **Dotmatics**, **LabArchives**. This is where experiments, structures and
  results are logged — and what Chemora's "remember your lab's history" story
  must integrate with.
- **Retrosynthesis / AI:** **Reaxys** and **SciFinder** retrosynthesis,
  **IBM RXN**, **ASKCOS** (MIT, open-source).
- **Analytical data:** **Mnova** (Mestrelab) for NMR/MS interpretation.

**Key insight for the India/mid-tier angle:** large Indian pharma and top CROs can
afford ChemDraw + Reaxys/SciFinder. The ~9,000 smaller companies the founder
described typically **cannot** — they often work in Word/PDF/paper and pirated or
shared single licenses. A genuinely affordable, locally-deployable tool that turns
their existing PDFs/notebooks into searchable, editable structures is a real,
underserved need. The on-prem/"your IP never leaves" angle also matters because
these firms are protective of proprietary routes.

---

## 5. UI/UX patterns worth borrowing

Study these before finalizing Chemora's UX:

- **DECIMER.ai** — the cleanest *free* "upload → segment → recognise → SMILES"
  flow. Good model for the extraction step's simplicity.
- **Reaxys** — gold standard for presenting an extracted structure alongside its
  data (conditions, yields, references). Borrow the "structure + evidence
  side-by-side" layout.
- **IBM RXN for Chemistry** — modern, friendly UX for a chemistry-AI web app;
  good reference for making something technical feel approachable.
- **ChemDraw / Ketcher canvas** — the editing surface itself (CanvasPrime already
  mirrors this).
- **The winning pattern for Chemora Extract:** *left = original PDF page,
  right = extracted structure on an editable canvas.* Source and result together,
  correction in place. That split **is** the product — it turns "trust the AI
  blindly" (scary for a chemist) into "AI did the work, I approve it" (easy yes).

---

## 6. Where Chemora's real opportunity is

Putting it together, the defensible positioning:

1. **Affordable + self-serve** where Reaxys/ChemDraw are enterprise-priced.
2. **PDF/notebook → editable canvas**, the workflow no single tool packages for
   small labs.
3. **On-prem / "your IP stays yours"**, which the target customers care about and
   the cloud giants can't easily match.
4. **Human-in-the-loop by design** — honest about accuracy, which builds the trust
   these firms need before adopting.
5. **Build on open-source OCSR** (OpenChemIE/DECIMER/MolScribe) + an open canvas
   (Ketcher) → low cost base → affordable price point is actually possible.

The prototype in this repo demonstrates #2 and #4 end-to-end, on an open-source
stack, which is the proof the rest of the positioning rests on.

---

## References (for the founder / for follow-up)

- OpenChemIE — github.com/CrystalEye42/OpenChemIE ; paper arXiv:2404.01462
- DECIMER.ai — decimer.ai ; DECIMER PyPI `decimer` ; Nat. Commun. 14:5045 (2023)
- MolScribe — github.com/thomas0809/MolScribe ; JCIM 63(7):1925 (2023)
- RxnScribe — github.com/thomas0809/RxnScribe
- RDKit — rdkit.org
- Ketcher — github.com/epam/ketcher
- Reaxys — elsevier.com/solutions/reaxys
- ChemDraw (Signals) — revvitysignals.com
- ChemAxon Document to Structure — docs.chemaxon.com
- IBM RXN — rxn.res.ibm.com
- ASKCOS (open retrosynthesis) — askcos.mit.edu
