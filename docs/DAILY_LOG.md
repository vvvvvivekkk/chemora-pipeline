# Daily log

One entry per working day, newest first. Numbers, not adjectives. The "learned"
line is written by the developer in his own words, not by Claude.

Template:

```
## YYYY-MM-DD
- Did:
- Score before → after:  (copy the scoreboard lines)
- New failures:          (→ FAILURE_LOG.md)
- Blockers:
- Learned (own words):
- Tomorrow:
```

---

## 2026-11-02
- Did: received the 20-step brief from SciChem AI; wrote `CLAUDE.md` (rules),
  `docs/PLAN.md` (phases), this log, `docs/FAILURE_LOG.md`. Cycle 7 probe
  (substitute text for R′/R″ picture-labels, 51 crops on CPU) was running on the
  local machine when it went offline — result not yet read.
- Score before → after: unchanged; last line is cycle 6, commit `4cc4f43`:
  structures exact 44/52 (84.6%), reactions exact 7/10, false positives 0 / 1,
  synthetic 35/47.
- New failures: none today (no model runs).
- Blockers: machine offline mid-probe; cycle 7 result pending.
- Learned (own words): _[Vivek fills this in — e.g. what "valid ≠ correct"
  means, why the answer key must not come from the model]_
- Tomorrow: read cycle 7 result, keep/revert, log it; start Phase 0 learning
  items (pipeline on paper, 5 SMILES ↔ pictures).

## 2026-10-31 (reconstructed from scoreboard and commit history)
- Did: cycles 1–6 of post-processing fixes on paper bjoc-20-228.
- Score: 26/52 → 36 (cycle 1, NMe₂ expansion) → 40 (cycle 2, picture-label
  R-groups → `*`) → 40 + reactions 2→4/10 (cycle 3, reaction parts through
  `postprocess_smiles`) → 37 (cycle 4, plain-R input substitution, **reverted**)
  → 42 (cycle 5, bond-order repair with ≤2-atom guard) → 44/52, reactions 7/10
  (cycle 6).
- New failures: see FAILURE_LOG.md entries F-001 … F-006.
- Blockers: none.
- Learned: fixing the model's *output* with code beat changing its *input*
  image (cycle 4 lesson).
