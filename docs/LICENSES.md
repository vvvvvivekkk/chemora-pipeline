# Licence table (Phase 1, brief step 9)

Rule: nothing marked **unsure** ships in the product until resolved. "Weights"
means the trained model files, which are often licensed separately from the code.

| Tool | Used for | Code licence | Weights licence | Commercial use | Verified on | Source |
|---|---|---|---|---|---|---|
| RDKit | validate / canonicalize / draw | BSD-3-Clause | — | yes | _todo_ | https://github.com/rdkit/rdkit |
| Streamlit | UI | Apache-2.0 | — | yes | _todo_ | https://github.com/streamlit/streamlit |
| Ketcher (via streamlit-ketcher) | editable canvas | Apache-2.0 | — | yes | _todo_ | https://github.com/epam/ketcher |
| OpenChemIE (MolDetect, MolScribe, RxnScribe code) | detection, reading, reactions | MIT | **unsure** — check each HF model card | **unsure** | _todo_ | https://github.com/CrystalEye42/OpenChemIE |
| MolScribe weights | structure reading | — | **unsure** | **unsure** | _todo_ | https://huggingface.co/yujieq/MolScribe |
| RxnScribe weights | reaction schemes | — | **unsure** | **unsure** | _todo_ | https://huggingface.co/yujieq/RxnScribe |
| MolDetect weights | structure detection | — | **unsure** | **unsure** | _todo_ | (inside OpenChemIE) |
| PyMuPDF (`fitz`) | PDF → page images | **AGPL-3.0 or paid commercial** | — | **no, unless licence bought** | _todo_ | https://pymupdf.readthedocs.io/en/latest/about.html |
| pypdfium2 (candidate replacement) | PDF → page images | Apache-2.0 / BSD-3 | — | yes | _todo_ | https://github.com/pypdfium2-team/pypdfium2 |
| poppler / pdftotext | text layer, HF Spaces system dep | GPL-2/3 | — | ok as a separate process; check | _todo_ | https://poppler.freedesktop.org |
| PyTorch 1.13 | model runtime | BSD-3 | — | yes | _todo_ | https://github.com/pytorch/pytorch |
| DECIMER (Phase 3) | second structure reader | MIT | **unsure** | **unsure** | _todo_ | https://github.com/Kohulan/DECIMER-Image_Transformer |
| MolNexTR (Phase 3) | third structure reader | **unsure** | **unsure** | **unsure** | _todo_ | https://github.com/CYF2000127/MolNexTR |
| Tesseract (Phase 4) | conditions OCR | Apache-2.0 | Apache-2.0 | yes | _todo_ | https://github.com/tesseract-ocr/tesseract |
| PaddleOCR / EasyOCR (Phase 4) | conditions OCR, second engine | Apache-2.0 | **check** | check | _todo_ | — |
| OPSIN | building answer keys only (not shipped) | Artistic-2.0 | — | not in product | _todo_ | https://opsin.ch.cam.ac.uk |

Action for the developer: open each link, read the LICENSE file and the model
card, fill "Verified on" with the date, and change every **unsure**.
The PyMuPDF row is the one most likely to force a code change.
