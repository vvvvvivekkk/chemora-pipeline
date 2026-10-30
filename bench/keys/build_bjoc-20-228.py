"""Build the answer key for bjoc-20-228 (Beilstein J. Org. Chem. 2024, 20, 2708).

Every SMILES comes from the paper, never from the model:
  * named compounds -> OPSIN name-to-structure on the paper's own IUPAC names
    (experimental section) — `opsin` below, fetched once and pasted here;
  * unnamed / generic drawings -> read from the scheme, using a substitution
    template that is first checked against every OPSIN-named member of its
    series, so a template error cannot slip in silently.
R, R1, R2, R', R'', X in a drawing -> '*'.
Boxes are normalised (x1, y1, x2, y2) page regions, read from the page images.

    python bench/keys/build_bjoc-20-228.py   # -> bench/keys/bjoc-20-228.json
"""
import json
from pathlib import Path

from rdkit import Chem

# OPSIN output for the names in the paper (opsin.ch.cam.ac.uk, fetched 2026-10-31)
opsin = {
    "1": "ClC1=CC(N(C(N1C)=O)C)=O",                 # 6-chloro-1,3-dimethylpyrimidine-2,4(1H,3H)-dione
    "2": "BrC=1C(N(C(N(C1Cl)C)=O)C)=O",             # 5-bromo-6-chloro-1,3-dimethyl... (Scheme 1, ref 65)
    "4g": "CN1C(N(C(C(=C1C#CC1=CC=CC=C1)C1=CC=CC=C1)=O)C)=O",
    "4h": "FC1=CC=C(C=C1)C=1C(N(C(N(C1C#CC1=CC=C(C=C1)F)C)=O)C)=O",
    "4i": "S1C2=C(C(=C1)C=1C(N(C(N(C1C#CC1=CC=C(C=C1)C)C)=O)C)=O)C=CC=C2",
    "5a": "CN1C(N(C=2C=C(C3=C(C2C1=O)C=CC(=C3)C)C3=CC=C(C=C3)C)C)=O",
    "5d": "CN(C)C1=CC=C(C=C1)C=1C2=C(C=3C(N(C(N(C3C1)C)=O)C)=O)C=CC(=C2)C",
    "5f": "CN(C)C1=CC=C(C=C1)C=1C2=C(C=3C(N(C(N(C3C1)C)=O)C)=O)C=CC(=C2)C(F)(F)F",
    "5g": "CN1C(N(C=2C=C(C3=C(C2C1=O)C=CC=C3)C3=CC=CC=C3)C)=O",
    "5h": "FC1=CC=C(C=C1)C=1C2=C(C=3C(N(C(N(C3C1)C)=O)C)=O)C=CC(=C2)O",
    "5i": "CN1C(N(C=2C=C(C3=C(C2C1=O)C1=C(S3)C=CC=C1)C1=CC=C(C=C1)C)C)=O",
}


def canon(smi):
    m = Chem.MolFromSmiles(smi)
    assert m is not None, smi
    return Chem.MolToSmiles(m)


def sub(template, **groups):
    """Fill [A]/[B] slots; an empty group means H (the slot and its bond vanish)."""
    s = template
    for k, v in groups.items():
        s = s.replace(f"[{k}]", f"({v})" if v else "")
    return canon(s)


# 1,3-dimethyl-6-(4-A-phenylethynyl)-5-(4-B-phenyl)uracil  (series 4)
T4 = "c1[A]ccc(C#Cc2n(C)c(=O)n(C)c(=O)c2-c2ccc[B]cc2)cc1"
# 2,4-dimethyl-6-(4-A-phenyl)-8-B-benzo[f]quinazoline-1,3-dione  (series 5)
T5 = "c1[A]ccc(-c2cc3n(C)c(=O)n(C)c(=O)c3c3ccc[B]cc23)cc1"
assert sub(T4, A="", B="") == canon(opsin["4g"])
assert sub(T4, A="F", B="F") == canon(opsin["4h"])
assert sub(T5, A="C", B="C") == canon(opsin["5a"])
assert sub(T5, A="N(C)C", B="C") == canon(opsin["5d"])
assert sub(T5, A="N(C)C", B="C(F)(F)F") == canon(opsin["5f"])
assert sub(T5, A="", B="") == canon(opsin["5g"])
assert sub(T5, A="F", B="O") == canon(opsin["5h"])

ME, NME2, CF3 = "C", "N(C)C", "C(F)(F)F"
c4 = {k: sub(T4, A=a, B=b) for k, (a, b) in {
    "4a": (ME, ME), "4b": (ME, NME2), "4c": (ME, CF3),
    "4d": (NME2, ME), "4e": (NME2, NME2), "4f": (NME2, CF3)}.items()}
c4.update({k: canon(opsin[k]) for k in ("4g", "4h", "4i")})
c5 = {k: sub(T5, A=a, B=b) for k, (a, b) in {
    "5b": (ME, NME2), "5c": (ME, CF3), "5e": (NME2, NME2)}.items()}
c5.update({k: canon(opsin[k]) for k in ("5a", "5d", "5f", "5g", "5h", "5i")})

# generic drawings (R groups -> *), read from the schemes
G = {
    "4": sub(T4, A="*", B="*"),                                    # Scheme 1-3, Fig. 1 "previous work"
    "5": sub(T5, A="*", B="*"),                                    # Scheme 1, 3, Fig. 1 "this work"
    "3": canon("*c1ccc(C#Cc2n(C)c(=O)n(C)c(=O)c2Br)cc1"),
    "alkyne": canon("*c1ccc(C#C)cc1"),                             # R'-phenylacetylene
    "boronic": canon("*c1ccc(B(O)O)cc1"),                          # R''-phenylboronic acid
    "F1_6aryl": canon("*c1ccc(-c2cc(=O)n(*)c(=O)n2*)cc1"),          # Miyasaka 1986 [62]
    "F1_5aryl": canon("*c1ccc(-c2cn(*)c(=O)n(*)c2=O)cc1"),          # Bristol 1905 [64]
    "F1_6alkynyl": canon("*c1ccc(C#Cc2cc(=O)n(C)c(=O)n2C)cc1"),     # Miyasaka 1986 [63]
    "F1_5alkynyl": canon("*c1ccc(C#Cc2cn(*)c(=O)n(*)c2=O)cc1"),     # Barr 1981 [51]
    "A": canon("O=C1N(*)C(=O)N(*)C2=C1c1ccccc1*2"),                # Fig. 2 A (X -> *)
    "B": canon("O=c1n(*)c(=O)c2cc3ccccc3nc2n1*"),                  # Fig. 2 B
    "C": canon("*n1c(=O)c2c3ccc(*)cc3c(-c3ccc(*)cc3)cc2n(*)c1=O"),  # Fig. 2 C
    # Scheme 4 intermediates (vinyl cations, protonated anilines)
    "Ia": canon("C[NH+](C)c1ccc([C+]=Cc2n(C)c(=O)n(C)c(=O)c2-c2ccc(C)cc2)cc1"),
    "Ib": canon("Cc1ccc([C+]=Cc2n(C)c(=O)n(C)c(=O)c2-c2ccc([NH+](C)C)cc2)cc1"),
}
assert G["C"] == sub(T5.replace("n(C)", "n(*)"), A="*", B="*")

S = []  # structures: (page, label, smiles, box)


def add(page, label, smi, box):
    S.append({"id": f"p{page}-{label}", "page": page, "label": label,
              "smiles": smi, "box": box, "junk": False})


def junk(page, label, box):
    S.append({"id": f"p{page}-{label}", "page": page, "label": label,
              "smiles": None, "box": box, "junk": True})


junk(1, "journal-logo", [0.09, 0.05, 0.18, 0.115])
junk(1, "cc-by-badge", [0.49, 0.46, 0.60, 0.50])
# page 2 — Figure 1 and Figure 2
add(2, "F1-Miyasaka62", G["F1_6aryl"], [0.11, 0.12, 0.25, 0.25])
add(2, "F1-Bristol64", G["F1_5aryl"], [0.26, 0.15, 0.46, 0.25])
add(2, "F1-Miyasaka63", G["F1_6alkynyl"], [0.11, 0.27, 0.25, 0.42])
add(2, "F1-Barr51", G["F1_5alkynyl"], [0.26, 0.31, 0.47, 0.42])
add(2, "F1-previous-work(4)", G["4"], [0.49, 0.22, 0.68, 0.38])
add(2, "F1-this-work(5)", G["5"], [0.69, 0.23, 0.89, 0.38])
add(2, "A", G["A"], [0.17, 0.50, 0.36, 0.61])
add(2, "B", G["B"], [0.36, 0.50, 0.56, 0.61])
add(2, "C", G["C"], [0.57, 0.50, 0.82, 0.61])
# page 3 — Scheme 1
add(3, "1", canon(opsin["1"]), [0.10, 0.16, 0.21, 0.29])
add(3, "2", canon(opsin["2"]), [0.30, 0.16, 0.42, 0.29])
add(3, "R'-alkyne", G["alkyne"], [0.42, 0.13, 0.49, 0.25])
add(3, "3", G["3"], [0.52, 0.15, 0.64, 0.32])
add(3, "R''-boronic", G["boronic"], [0.60, 0.13, 0.68, 0.24])
add(3, "4", G["4"], [0.71, 0.15, 0.90, 0.32])
add(3, "5", G["5"], [0.41, 0.33, 0.60, 0.50])
# page 4 — Scheme 2
add(4, "3", G["3"], [0.27, 0.16, 0.40, 0.32])
add(4, "R''-boronic", G["boronic"], [0.40, 0.14, 0.48, 0.25])
add(4, "4", G["4"], [0.52, 0.16, 0.70, 0.32])
for row, labels in enumerate([("4a", "4b", "4c"), ("4d", "4e", "4f"), ("4g", "4h", "4i")]):
    y1 = [0.34, 0.51, 0.69][row]
    for col, lab in enumerate(labels):
        x1 = [0.18, 0.40, 0.62][col]
        add(4, lab, c4[lab], [x1, y1, x1 + 0.20, y1 + 0.16])
# page 6 — Scheme 3
add(6, "4", G["4"], [0.19, 0.11, 0.40, 0.30])
add(6, "5", G["5"], [0.57, 0.12, 0.78, 0.30])
for row, labels in enumerate([("5a", "5b", "5c"), ("5d", "5e", "5f"), ("5g", "5h", "5i")]):
    y1 = [0.30, 0.48, 0.64][row]
    for col, lab in enumerate(labels):
        x1 = [0.17, 0.39, 0.61][col]
        add(6, lab, c5[lab], [x1, y1, x1 + 0.21, y1 + 0.16])
# page 7 — Scheme 4 and Figure 3
add(7, "4d", c4["4d"], [0.14, 0.11, 0.32, 0.28])
add(7, "Ia", G["Ia"], [0.41, 0.11, 0.64, 0.30])
add(7, "5d(4a)", c5["5d"], [0.72, 0.11, 0.90, 0.28])
add(7, "4b", c4["4b"], [0.14, 0.29, 0.34, 0.44])
add(7, "Ib", G["Ib"], [0.41, 0.29, 0.64, 0.45])
for lab, x1, x2 in [("5a", 0.09, 0.24), ("5d", 0.23, 0.37), ("5f", 0.36, 0.50),
                    ("5g", 0.49, 0.61), ("5h", 0.60, 0.75), ("5i", 0.74, 0.90)]:
    add(7, f"Fig3-{lab}", c5[lab], [x1, 0.49, x2, 0.65])
junk(7, "spectra-plots", [0.13, 0.66, 0.86, 0.88])

R = [  # reactions: reactants >> products, reagents drawn as structures count as reactants
    {"id": "S1-i", "page": 3, "reactants": [canon(opsin["1"])], "products": [canon(opsin["2"])]},
    {"id": "S1-ii", "page": 3, "reactants": [canon(opsin["2"]), G["alkyne"]], "products": [G["3"]]},
    {"id": "S1-iii", "page": 3, "reactants": [G["3"], G["boronic"]], "products": [G["4"]]},
    {"id": "S1-iv", "page": 3, "reactants": [G["4"]], "products": [G["5"]]},
    {"id": "S2", "page": 4, "reactants": [G["3"], G["boronic"]], "products": [G["4"]]},
    {"id": "S3", "page": 6, "reactants": [G["4"]], "products": [G["5"]]},
    {"id": "S4a-1", "page": 7, "reactants": [c4["4d"]], "products": [G["Ia"]]},
    {"id": "S4a-2", "page": 7, "reactants": [G["Ia"]], "products": [c5["5d"]]},
    {"id": "S4b", "page": 7, "reactants": [c4["4b"]], "products": [G["Ib"]]},
]

out = Path(__file__).with_name("bjoc-20-228.json")
out.write_text(json.dumps({
    "paper": "bjoc-20-228",
    "citation": "de Abreu, Ehlers, Langer, Beilstein J. Org. Chem. 2024, 20, 2708-2719, doi:10.3762/bjoc.20.228",
    "structures": S, "reactions": R}, indent=1))
print(f"{out.name}: {sum(not s['junk'] for s in S)} structures, "
      f"{sum(s['junk'] for s in S)} junk regions, {len(R)} reactions")
