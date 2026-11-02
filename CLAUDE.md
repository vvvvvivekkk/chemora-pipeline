# Chemora Extract — rules for every Claude session

Read this before doing anything in this repo. These rules come from the project
owner (SciChem AI) and override convenience.

## Hard rules (never break)

1. **This is an offline pipeline. No LLM may generate SMILES.** Structures come
   only from dedicated recognition models (MolScribe, RxnScribe, DECIMER,
   MolNexTR). Post-processing may be ordinary code (string/graph edits, RDKit),
   never a language model guessing a molecule. If a step would need an LLM to
   decide what a molecule is, stop and say so.
2. **Isolated environments, pinned versions.** Never install into the global
   environment. Each model family gets its own venv/conda env with every version
   pinned. Model weights are never committed (see `.gitignore`).
3. **Verify before claiming.** Nothing "works" until `python bench/evaluate.py
   bjoc-20-228` (and the synthetic benchmark it runs) has been re-scored and the
   new line is on `bench/scoreboard.md`. Report the numbers, not the feeling.
4. **Fail loudly.** If something errors, stop and show the exact error. Do not
   work around it silently, do not weaken a test or a guard to make it pass, do
   not skip a failing case.
5. **Explain chemistry terms.** The developer is a chemistry beginner. Define
   every new term simply, with one example, the first time it appears.

## Working rules

- **Ask before changing core logic:** voting/agreement, validation, the
  confidence rules, the benchmark scoring, or the answer keys.
- **Branch per change, small commits**, each message saying what changed and
  why. `main` only receives work that is scored and not worse than the previous
  scoreboard line.
- **One cycle = one failure category.** Probe first (measure on cached outputs),
  then fix, then score. Keep only if the target improves and nothing else
  regresses on the same score line. A reverted attempt stays in git history and
  gets a `docs/FAILURE_LOG.md` entry.
- **RDKit "valid" ≠ correct.** Valid only means parseable. Never label ✅ as
  "correct" in code, UI text, or reports.
- **Answer keys** come from verified sources only (OPSIN on the paper's own
  IUPAC names, PubChem, ChemDraw export) — never from the model, never typed
  from memory, never from an LLM.
- **Logs are kept in the repo:** `docs/DAILY_LOG.md` (daily, in numbers) and
  `docs/FAILURE_LOG.md` (every wrong read worth remembering). Update both at the
  end of each session.
- The plan is `docs/PLAN.md`. Work the phases in order unless the owner says
  otherwise.

## Scoring commands

```bash
python bench/evaluate.py bjoc-20-228      # paper + synthetic, appends scoreboard line
python bench/benchmark.py run             # synthetic only
python test_e2e.py                        # smoke test
```
