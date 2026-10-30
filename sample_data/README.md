# Sample test PDFs

You need a few chemistry PDFs with clear structure/reaction figures to test with.
Don't commit large PDFs to the repo — just grab a couple locally for the demo.

## Good sources (open access, legal to download)

- **Open Reaction Database (ORD)** — https://open-reaction-database.org
  Real reaction data; good for reaction-scheme testing.
- **PubChem** — https://pubchem.ncbi.nlm.nih.gov
  Single-molecule pages; export/print a structure to PDF for a clean easy test.
- **ChemRxiv** — https://chemrxiv.org
  Open-access chemistry preprints with real figures.
- **PLOS / MDPI / Beilstein J. Org. Chem.** — open-access journals full of
  reaction schemes.

## Tips for a convincing demo

1. **Start easy.** A page with ONE clear, printed structure recognises best.
   Prove the pipeline works, then show a harder multi-structure page.
2. **Printed > hand-drawn.** These models handle clean printed figures far better
   than hand-drawn ones.
3. **Pre-run before the demo** so model weights are cached and it's fast live.
4. Keep a **known-good SMILES** handy (e.g. aspirin `CC(=O)Oc1ccccc1C(=O)O`) to
   show the `manual` engine + canvas instantly if the network is slow.
