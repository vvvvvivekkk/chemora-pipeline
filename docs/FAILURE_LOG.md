# Failure log

Every misread or failed attempt worth remembering. Format per the brief:
image, model answer, correct answer, probable reason. Plus: status (open /
fixed in cycle N / won't fix) so the same lesson is not relearned.

Paper = Beilstein J. Org. Chem. 2024, 20, 2708 (bjoc-20-228) unless stated.
Key: `bench/keys/bjoc-20-228.json`. Report with crops: `bench/report.html`
(regenerate with `python bench/evaluate.py bjoc-20-228`).

| id | image / crop | model answer | correct answer | probable reason | status |
|---|---|---|---|---|---|
| F-001 | NMe₂ substituent on 5d, 5f | label read as R-group `*` | `N(C)C` | MolScribe treats unexpanded text labels at a ring edge as generic R | fixed, cycle 1 (label expansion) |
| F-002 | R′ circle / R″ box picture-labels (series 4, 5) | `Pd`, `Rn`, `Pa`, `Rf`, `Fl`, `[2H]` atoms | `*` | labels are pasted images, not text; model guesses a rare element | fixed, cycle 2 (picture-label → `*` by position) |
| F-003 | reactant/product crops inside reactions (S1-ii, S3, S1-iii) | same raw misreads as F-002 | post-processed SMILES | RxnScribe calls MolScribe itself; results bypassed `postprocess_smiles` | fixed, cycle 3 |
| F-004 | generic 4 on p2, p3, p6 (after plain-"R" input substitution) | `P` (phosphorus), sometimes `C` or nothing | `*` | a bare "R" glyph next to a bond stub looks like "P" to the model; changing input image swapped one misread for another | attempt reverted, cycle 4 — **lesson: fix output, not input** |
| F-005 | 5i (p6, p7 Figure 3), fused benzothiophene | RDKit invalid: two carbons with valence 5 | `CN1C(N(C=2C=C(C3=C(C2C1=O)C1=C(S3)C=CC=C1)C1=CC=C(C=C1)C)C)=O` | atoms and connectivity right; double bonds placed wrongly in the fused ring | fixed, cycle 5 (re-derive bond orders; guard: ≤2 over-valent atoms, each +1) |
| F-006 | generic 5 (p2, p6 invalid; p3 near) | 51-atom garbled graph (26 expected), 8–9 over-valent atoms | series-5 template with `*` | two picture-labels *inside* a fused ring system confuse the reader; too broken for the cycle-5 repair guard (correctly) | open — cycle 7 probing substitute labels on cached crops |
| F-007 | Bristol, Barr (benzene rings) | 7-membered ring | 6-membered benzene | unknown; drawing style (coloured highlight?) — probe needed | open |
| F-008 | 4h, 5h: F at top of ring | `E` → becomes `*` | `F` | "E" is in MolScribe's R-group symbol list, so a misread F silently turns into an R-group | open — candidate rule: `*` at a position where the key says F is suspicious; needs a general test, not a key-specific hack |
| F-009 | Ia, Ib vinyl cations; reactions S4a-1, S4a-2, S4b | `+` charge missed | cation | the `+` is a separate small glyph beside the carbon; unclear if the model sees it at all | partly addressed cycle 6 (reactions 4→7/10); check which remain |
| F-010 | X scaffold "A" | wrong | — | not yet investigated | open |
| F-011 | synthetic paclitaxel (11 stereocentres) | over-valent by 2 | — | very large molecule; reader garbles; correctly left flagged by cycle-5 guard | won't fix via post-processing; voting/second reader may help |
| F-012 | bjoc.20.233 CF₂ drawn as stacked C with F₂ above/below | `CH₂CH₂Br` (fluorines dropped) | `CF₂CF₂Br` | stacked labels lose the stacked part | open — needs its own key (Phase 2) |
| F-013 | bjoc.20.233 printed `*` marking a new stereocentre | read as R-group | plain carbon | `*` glyph is literally the R-group symbol | open |
| F-014 | journal logo, licence badge (before cleanup) | `Co.Co.Co.Co` etc. | nothing (not a molecule) | MolDetect boxes non-chemistry graphics | fixed (cleanup drops carbon-free output) |
