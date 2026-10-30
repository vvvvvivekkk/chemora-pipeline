"""End-to-end check: sample PDF -> pipeline -> validated SMILES.

Run inside the `chemora` env:  python test_e2e.py
Fails (non-zero exit) if the pipeline errors or finds no valid molecule.
"""
import sys
import time

from pipeline import config
from pipeline.orchestrator import process_pdf


def main() -> int:
    pdf = open("sample_data/sample_molecules.pdf", "rb").read()
    t0 = time.time()
    result = process_pdf(pdf, progress_callback=lambda f, m: print(f"[{f:4.0%}] {m}", flush=True))
    dt = time.time() - t0

    print(f"\nengine={result.engine}  pages={len(result.pages)}  "
          f"molecules={result.total_molecules}  reactions={result.total_reactions}  "
          f"({dt:.0f}s)")
    for p in result.pages:
        for m in p.molecules:
            print(f"  page {p.index + 1}: valid={m.is_valid}  smiles={m.canonical_smiles or m.raw_smiles!r}  {m.message}")
        for r in p.reactions:
            print(f"  page {p.index + 1}: reaction {r}")

    assert result.pages, "no pages rendered"
    if config.ENGINE == "manual":
        print("manual mode: pages rendered, no model run — OK")
        return 0
    valid = [m for p in result.pages for m in p.molecules if m.is_valid]
    assert valid, "no valid molecules extracted"
    print(f"\nE2E OK: {len(valid)} valid molecule(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
