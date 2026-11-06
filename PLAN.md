# Agentic solutions to the MDS CMSDAS exercises (revision 2)

## Context
Repeat the earlier plan (`i-am-trying-to-lively-snail.md`) with three changes requested by the user:
1. The first result filled code blanks but did not answer the **discussion questions** in the episodes. Every `Discussion N.M`, `Question N.M`, `Exercise`, and `Check your understanding` block must now get an explicit written answer.
2. Agents may read the **`gh-pages` branch** (episode text) and nothing else outside the allowed sources below.
3. **Never touch** `/uscms_data/d3/kkwok/cmsdas/MDS_Ex/CMSSW_14_1_0_pre4/src/solution` (no read, write, ls, grep, copy). This is different from the repo-local `solutions_agentic/`. Agents must be told this explicitly in their prompts.

Existing state: branch `agentic-solution-2026-11-06-opencode-qwen3` (untracked `.opencode/`) and `agentic-solution-2026-11-06-claude-code-haiku-5-5`. Both are prior attempts; leave them untouched. `solutions_agentic/ep02..ep07*.md` exist on the current branch from the first result.

## Allowed sources (everything else off-limits)
- `gh-pages` branch, read-only via `git show gh-pages:_episodes/NN-*.md` (episodes 02-07; 01 for background). Do not check it out.
- `2025` skeleton branch / current working tree: `scripts/*.ipynb`, `scripts/run_combine.py`, `lib/{helper,util,tdrstyle}.py`.
- Own physics knowledge.
- **Forbidden:** `master`, the `solution` directory above, `.opencode/`, other agentic branches, any other path on disk.

## Branch and layout
- New branch: `git switch -c agentic-solution-<date>-<tool>-<model> 2025` (e.g. `...-claude-code-haiku-5-5`; pick a name not already used; no merge/diff against master). This plan is model-agnostic: whichever model/tool executes it fills in its own tag and reports its own token usage.
- Deliverables:
  1. `solutions_agentic/epNN_*.md` per episode, restructured around the episode's own discussion flow: for each exercise, in order, (a) filled code, (b) student-level explanation, (c) **an answer to each Discussion/Question/Exercise block quoted verbatim from gh-pages**, with reasoning grounded in the notebook plots/numbers and physics (e.g. lifetime reweighting, geometric acceptance, cluster variables, ABCD independence/closure, limit scaling with ctau). Plus `README.md` with a coverage checklist and token estimate.
  2. In-place fills of `scripts/*.ipynb` and `scripts/run_combine.py`, with a markdown cell under each solution containing the explanation and the discussion answers. Edit JSON via script; nothing executed; outputs left as is.
- Copy this plan verbatim into the repo directory as `PLAN.md` (the "this directory" for re-use) once approved. Confirm location with the user if they meant another directory.

## Process
1. Pull the question inventory: for each episode run `git show gh-pages:_episodes/NN-*.md` and extract every Discussion/Question/Exercise/Check block with its line numbers (known: Disc 2.1-2.3; 3.1-3.3 plus the two cluster-property exercises; 4.1-4.5, Question 4.1, 4.7; 5.1; 6.1, Question 6.1, 6.1 (limits), final Exercise; ep07 questions). Write the inventory to the scratchpad as the coverage checklist.
2. Extract notebooks to text with outputs stripped (ABCD_validation is 225 KB, mostly output; strip first). Output text and plots are still useful for answers, so keep a short summary of key printed numbers where they back a discussion answer.
3. Up to 3 parallel agents ({02,03}, {04,05}, {06,07}), each given: the episode text, stripped notebooks, lib helpers, the checklist slice, and the forbidden-path rule. Each returns code fills and an answer for every checklist item. Main agent reviews for consistency (variable and branch names, cm/m units, cuts quoted in the lesson text).
4. Apply fills via a JSON-editing script; write markdown; syntax check only (`json.load`, `python -m py_compile`).
5. Commit on the new branch (with whatever attribution the running tool requires); no push.

## Token measurement (for cross-model comparison)
- Rough content volume, independent of model: episodes ~12k tokens to read, notebooks and lib ~25-40k, written output ~45-65k, subagent overhead ~60-100k. Tokenizers differ, so treat as order of magnitude.
- The measured figure is the executing tool's own usage report (e.g. `/cost` in Claude Code), recorded in `README.md` together with: model name, tool, number of subagents, wall-clock time, and the date. Use the same task scope and checklist for every model so runs are comparable.

## Verification
- Coverage: every checklist item has an answer in the md and in the notebook; README lists them.
- Notebooks parse as nbformat JSON; `run_combine.py` compiles.
- `git log`/diff shows branch based on 2025, no master or `solution`-directory content, and no files outside `scripts/`, `solutions_agentic/`, `PLAN.md` changed.
- Spot-check physics: decay-length reweighting formula, geometric acceptance, ABCD closure definition, Combine flags match episode 6 text.
