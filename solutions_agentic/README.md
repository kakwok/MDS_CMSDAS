# Agentic solutions: run record and coverage checklist

Branch: `agentic-solution-2026-11-06-claude-code-haiku-5-5-v2` (based on `2025`).
Plan: `PLAN.md` at the repo root (model-agnostic; re-run with any model).
Nothing was executed: no notebook, script, or Combine run. Checks were syntax only (`json.load`, `ast.parse`, `py_compile`).

## Run record (fill in for each model)

| Field | Value |
|---|---|
| Model | Haiku 5.5 (`claude-haiku-5-5`) |
| Tool | Claude Code |
| Subagents | 3 general-purpose workers, one per episode group ({02,03}, {04,05}, {06,07}) |
| Subagent tokens (reported) | 121,702 + 106,594 + 107,638 = 335,934 |
| Subagent tool uses | 18 + 16 + 14 = 48 |
| Main-session tokens | to be taken from `/cost` |
| Wall-clock (subagents) | 263 s + 207 s + 166 s (run in parallel) |
| Date | 2026-11-06 |

## Coverage checklist

Status key: **done** = answered in md and notebook; **partial** = answered, but the answer needs a run to confirm a number.

### Episode 02 (LLP lifetimes)
- Discussion 2.1, 2.2, 2.3: done (`ep02_LLP_lifetimes.md`)
- Ex0 to Ex3 code cells: done (LLP_lifetimes cells 7, 10, 13, 16, 19, 22; acceptance loop filled)
- Notebook question "Do you understand the shape?": done

### Episode 03 (MDS reconstruction)
- Discussion 3.1 (both sub-questions), 3.2, 3.3 (three parts): done (`ep03_MDS_reco.md`); 3.3 station-level answers are qualitative
- Challenges "compute cluster properties" and "plot cluster properties": done
- Notebook questions (cluster shape, fill-variable, inclusive efficiency): done
- Inclusive efficiency and station efficiency values: partial (need a run)

### Episode 04 (analysis strategy)
- Check your understanding, Discussion 4.1 to 4.5, 4.7, Question 4.1 with its solution: done (`ep04_analysis_strategy.md`)
- Notebook question cells (10): done
- Blank code cells (17): filled
- Trigger efficiency, veto and cut efficiencies: partial (need a run)

### Episode 05 (ABCD)
- Discussion 5.1, notebook correlation question, OOT checklist: done (`ep05_ABCD.md`)
- Blank code cells (4): filled; cell 5 `REPLACE_ME` set to `-12.5`
- Pearson rho value and ABCD agreement: partial (need a run)

### Episode 06 (limits)
- Discussion 6.1 (both blocks), Question 6.1 with its solution, 15/55 GeV exercise, Combine/limit section: done (`ep06_Limits.md`)
- `create_datacard` cells (ggH QCD scale, data and signal bin layout): filled
- `limitPlot` cell (paths and isfile guard; BR = r/sigA): filled
- `scripts/run_combine.py` (real card path, `--freezeParameters norm`, `--dryRun`): filled
- "Does background agree?" answer: partial (needs a Combine run)

### Episode 07 (event display)
- Checklist and shell/CernBox/cmsShow steps: done (`ep07_event_display.md`)
- `event_display` cell 7: fixed indexing (mask applied before indexing, `min(10, ...)`)

## Open issues (for review before merging)

1. **Unverified field names.** `rechits.Ring` (ME11/12 selection) and `signal.cscRechitClusterZ` (cell 69) are not in the lesson text, the lib helpers, or any saved notebook output. Check `rechits.fields` (MDS_reconstruction cell 15) and `signal.fields` (cell 48) on a real run.
2. **Event display threshold.** Code uses `Nhits > 500`; episode 07 refers to the signal region of 130 hits / 0.75. Kept 500; needs a decision.
3. **ABCD convention.** `ABCD_validation` cell 18 says `A = B x C / D`. The lesson model (and the datacard formula) gives `A = B x D / C`. The code and notes use `B x D / C`. The notebook's own text should be corrected.
4. **OOT cut wording.** The `ABCD_validation` cell text says `< 12.5 ns` but the instruction says to set `REPLACE_ME` to `-12.5` (the early out-of-time region). A note was added as a markdown cell; the text itself was not changed.
5. **Sentinel values.** The `-1` sentinel for empty events falls in the underflow bin and inflates efficiencies. The `|eta|` variable uses `-999` instead. Other variables that use `-1` should be checked.
6. **MC time window.** The worker used `-5 < t < 12.5 ns` for MC. This was a choice, not from the lesson.
7. **Endcap cuts (Ex3).** `|eta| < 2.4`, `568 < |z| < 1075`, `|r| < 700` were chosen from `lib/util.py` geometry. Acceptance (target 15 to 20%) was not verified.
8. **Guard added.** `computeStationProp` returns `-999` when no station has at least 10 hits. This guard is not in the template.
9. **Numerical answers.** No discussion answer quotes a number. Each one says what to look for.
10. **Leftover banner comment.** LLP_lifetimes Ex3 cell keeps its `FILL IN HERE` banner comment above the filled code. The code is complete.
