# Episode 06: Results and Statistical Analysis (Limits)

Scope: `scripts/create_datacard.ipynb`, `scripts/limitPlot.ipynb`, `scripts/run_combine.py`. The Higgs Combine steps are described but were not executed here (syntax checks only).

## Inventory of blocks (episode 06-Limits.md)

| # | Block | Line range | Status |
|---|---|---|---|
| 1 | Checklist: open `create_datacard.ipynb` | 23-26 | Covered by notebook fills below |
| 2 | Discussion 6.1 (datacards / ABCD constraint) | 29-33 | Answered below |
| 3 | Question 6.1 (challenge: `norm` rateParam) | 35-37 | Answered below |
| 4 | Solution 6.1 (provided) | 39-44 | Checked, consistent with the code |
| 5 | Results: background-only fit, "Does the background prediction agree with the observation?" | 51-87 | Answered below (needs a Combine run to give a number) |
| 6 | Limits: AsymptoticLimits command, `higgsCombineTest.AsymptoticLimits.mH120.root` | 91-104 | Answered below |
| 7 | Checklist: open `scripts/run_combine.py` | 106-109 | Filled (full file) |
| 8 | Checklist: open `limitPlot.ipynb` | 116-119 | Filled; asked to compare with slide 25 (not checkable here) |
| 9 | Discussion 6.1 (second, same label) (limit shape vs ctau) | 121-127 | Answered below |
| 10 | Exercise (15 and 55 GeV) | 129-133 | Answered below as a plan (not executed) |

Labels: "Discussion 6.1" occurs twice in the episode (lines 29 and 123). Both are answered separately.

Total: 10 items, of which 2 Discussion blocks, 1 Question/Challenge, 1 Solution, 1 Exercise, 3 Checklists, and 2 in-text questions.

## Code fills

### 1. `create_datacard.ipynb`, cell 8 (`ggH_qcdScale`)

```python
ggH_qcdScale[key] = [0.067] * 4 + [0.046]*4 #down*4/up*4
```

Explanation: the list has eight entries, the first four are the downward shifts for bins A, B, C, D and the last four are the upward shifts. The QCD scale effect is an overall normalization shift of the ggH cross section (6.7% down, 4.6% up), so the same value goes into every bin. `make_datacard` writes `1-down/1+up`, giving `0.933/1.046`, which matches the example datacard in the notebook.

Pitfall: the sign convention is down first, then up. Reversing them swaps the two systematic shifts.

### 2. `create_datacard.ipynb`, cell 12 (data bins B, C, D and signal bins C, D)

```python
b = np.count_nonzero(np.logical_and(Nhits[k] >= N_RECHIT_CUT, np.abs(deltaPhi[k]) >= DPHI_CUT))
c = np.count_nonzero(np.logical_and(Nhits[k] < N_RECHIT_CUT, np.abs(deltaPhi[k]) >= DPHI_CUT))
d = np.count_nonzero(np.logical_and(Nhits[k] < N_RECHIT_CUT, np.abs(deltaPhi[k]) < DPHI_CUT))
```

and, inside the signal loop:

```python
        # bin C
        cond = np.logical_and(Nhits[k]<N_RECHIT_CUT, deltaPhi[k]>=DPHI_CUT)
        signal_rate.append(np.sum(w[cond]))

        # bin D
        cond = np.logical_and(Nhits[k]<N_RECHIT_CUT, deltaPhi[k]<DPHI_CUT)
        signal_rate.append(np.sum(w[cond]))
```

Explanation: the ABCD plane is defined by two cuts, N_hits >= 130 and |dphi(cluster, MET)| < 0.75.
- A: high N_hits, low dphi (signal region). Already given.
- B: high N_hits, high dphi.
- C: low N_hits, high dphi (diagonal to A).
- D: low N_hits, low dphi.

The layout is chosen so that B and C share the high-dphi column, and C and D share the low-N_hits row. This gives `N_A = N_B * N_D / N_C` for the background, which is the relation written in the datacard (`(@0*@2/@1)`). Data and signal use the same four cuts, so the signal and observation are in the same bins.

Cross-check: in the example datacard, the signal yields are about 2010 (A), 7.9 (B), 0.84 (C), 286 (D). Most signal is in A and D (low dphi), and very little is in B and C (high dphi), which matches this layout.

Pitfalls:
- Use `>=` on N_hits in the signal region and `<` for the complement, so that every event falls in exactly one bin.
- The data counts use `np.abs(deltaPhi[k])`, as in the existing line for `a`. The signal bins use `deltaPhi[k]`, as in the existing bin B line. The two agree if dphi is non-negative.
- `ABCD_validation.ipynb` writes the relation as `A = B*C/D` with a different labelling. This is a different notebook and was not changed, but the datacard layout above follows the episode 05 text (N_A = c1*c2*Bkg_C) and the datacard formula. Check that the two notebooks use the same bin names before comparing numbers.

### 3. `limitPlot.ipynb`, cell 3 (paths and the limit line)

Changes:
- `dataCardDir = os.getcwd().replace('scripts', 'combine/datacards')`
- `limitDir = os.getcwd().replace('scripts', 'combine/limitTrees')`
- Added an `os.path.isfile(limitTrees[m])` guard before `uproot.open`, so a missing limit tree is skipped instead of raising.
- The line `limits[i] = np.array(T['limit'].array())/signal_norm` is kept as written. The comment that said "CHANGE THIS LINE" was replaced with a note on the normalization.

Explanation of the normalization (this is the point of the "CHANGE THIS LINE" comment):
- In the datacard, `signal_norm` is the signal yield in bin A at BR = 1 (the first `rate` entry).
- `run_combine.py` sets the Combine parameter `norm = 1/signal_norm`. Combine's expected signal is `signal_norm * norm * r`, so `r` equals the number of signal events in bin A.
- The BR limit is the event count divided by the BR = 1 yield: `BR_lim = r_lim / signal_norm`. This is `r_lim * norm_set`.

So the code's division by `signal_norm` is correct. Note that the wording "divide the limit by the norm you set" would give the wrong answer, because the norm you set is `1/signal_norm`. If you instead run Combine with `norm=1` (the manual command in the episode), `r` is already the BR and no division is needed.

`h_obs` (cell 5) is already present in the stripped notebook, so no change is needed there.

Pitfall: the `ctaus` list in `limitPlot.ipynb` has to match the card names written by `create_datacard.ipynb` and `run_combine.py`. The name pattern is `ggH_HToSSTobbbb_mh125_mx40_ctau{ct}mm`.

### 4. `scripts/run_combine.py` (full file)

Changes:
- `card_dir` now points to `<repo>/combine/datacards/`, derived from `__file__`. The placeholder `CMSSW_BASE + "path_to_datacard"` was replaced.
- `card_dir` has a trailing slash, so the file paths are `card_dir + name + '.txt'`.
- The Combine option is `--freezeParameters norm` (the episode uses this spelling). The old `--freezeParameter=norm` was changed.
- The `--dryRun` option, which was defined but unused, now prints the command and skips running it.
- Removed the unused `CMSSW_BASE` lookup.

The command it runs for each ctau is:

```
combine -M AsymptoticLimits <card_dir><name>.txt -n.<name> --setParameters norm=<1/signal_norm> --freezeParameters norm
```

Output: `higgsCombine.<name>.AsymptoticLimits.mH120.root`, moved into `combine/limitTrees/`. This is the file that `limitPlot.ipynb` reads from `limitDir`.

## Answers to the discussion and question blocks

### Discussion 6.1 (line 29-33, first occurrence)

> Do you understand the example datacards? Where is the ABCD relationship constrained for background?

Answer: yes, the four sections are the observation, the signal and background rates, the ABCD relations, and the signal uncertainties.
- The observation line gives the data counts in A, B, C, D.
- The background enters the `rate` line as a placeholder with value 1 (`bkg`). The actual background rate is set by the `rateParam` lines for B, C and D. Each of these is a free parameter in the range `[0, 7*observed]`. Their Poisson likelihood against the observed data is what constrains them.
- The A line `A rateParam chA bkg (@0*@2/@1)` sets the background in A to `B*D/C` (`@0`=B, `@1`=C, `@2`=D). This is the constraint: the background in A is not free, it is fixed by the other three bins.
- The signal uncertainties are log-normal (`lnN`) nuisances on the signal rate.

So the ABCD relationship is constrained in the `ABCD relationships` block, through the formula on the A line, with B, C, D floating against data.

### Question 6.1 (challenge, line 35-37) and Solution 6.1

> Do you know what does the rateParam `norm` do and why do we add the rateParam to change the normalization of the signal yield?

Answer: `norm rateParam * signal 1` multiplies the signal rate in every bin by a free parameter `norm`, starting at 1. We add it because the signal yield changes by about three orders of magnitude across lifetimes (O(1) to O(1000) events). Combine's signal strength `r` is only well-behaved when the fitted `r` is of order 1 to 15. So we scale the signal yield by `1/signal_norm` (through `norm`), keep `r` in that range, and multiply back at the end when we plot the limit.

The solution text is correct. One point to add: the `norm` is frozen at its set value (`--freezeParameters norm --setParameters norm=...`), so it is not fitted. Freezing is essential. If `norm` floated, it would be degenerate with `r`, and the limit on BR would be meaningless.

### Results question (line 87): "Does the background prediction agree with the observation?"

Not answered with a number. The background-only fit and the FitDiagnostics readout need a Combine run, which was not executed here.

What to check: run the two commands from the episode and compare `total_background` bin 1 (the A-bin prediction, from the fit with `r=0`) with `observation[0]` (the observed count in A). Agreement within the printed uncertainty means the ABCD prediction is consistent with data. The B, C and D bins are fitted to data, so only the A prediction is a real test.

### Discussion 6.1 (line 121-127, second occurrence)

> Did you expect the shape of the limit to look like this? Why does the sensitivity decreases for very long or very short lifetimes?

Answer: the shape is expected. The limit on BR is roughly the inverse of the number of signal events the analysis can see. The number of seen signal events is the lifetime times the probability that the LLP decays inside the muon system, which is at a distance of metres from the interaction point.
- Very short lifetimes (ctau of a few mm to cm): the LLP decays before the muon system, so few clusters are produced in the CSC. The signal yield is small, the BR limit is weak.
- Very long lifetimes (ctau of many metres): most LLPs escape the detector without decaying, again giving few signal clusters, and the limit is weak.
- Best sensitivity at ctau comparable to the detector distance, where the decay probability in the CSC is highest.

The plot in `limitPlot.ipynb` uses the log scale in ctau and BR, so the minimum of the expected curve should be near ctau of order 1 m for 40 GeV.

### Exercise (line 129-133): 15 and 55 GeV

> Now that you have the limit plot for 40 GeV, try to see if you can create datacards, run combine, and make limits plot for 15 and 55 GeV.
> Try to plot them in the same canvas, like the one in slide 25 of the introduction slide.

Plan (not executed):
1. In `create_datacard.ipynb`, set `mass = [15, 55]` (cell 2 and cell 8), and load the MS-15 and MS-55 signal ntuples, each with its own `OLD_CTAU` list. Also set the mass in the `modelName` string (`mx{m}`).
2. In `run_combine.py`, set `sample_name` per mass (for example `ggH_HToSSTobbbb_mh125_mx15_ctau...`). Loop over masses or run the script once per mass.
3. In `limitPlot.ipynb`, change the sample name `mx40` to the target mass in `limitTrees` and `dataCards`. Make one `TGraph` set per mass, and draw them on the same canvas with a legend entry per mass.
4. Check that the ctau list for each mass matches the signal samples available for that mass.

Pitfall: the ggH uncertainty dictionaries in `create_datacard.ipynb` are keyed by `MC_<mass>_<ctau>`, so they must be filled for every mass.

## Check status

- `python -I -m py_compile` passed on all new code cells (create_datacard cells 8 and 12, event_display cell 7, limitPlot cell 3) and on `scripts/run_combine.py`.
- The notebooks and the script were not executed. The Combine numbers and the slide 25 comparison were not checked.
