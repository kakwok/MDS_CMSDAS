# Episode 06: Results and limits (agentic solutions)

Scope: `scripts/create_datacard.ipynb`, `scripts/run_combine.py`, `scripts/limitPlot.ipynb`, and the exercise at the end of the episode.

## Create datacards (`scripts/create_datacard.ipynb`)

> Apply the signal region selections, add the systematic uncertainties, and write the signal and background yields in the ABCD plane to datacards.

### Blank 1: QCD scale uncertainty

```python
ggH_qcdScale[key] = [0.067] * 4 + [0.046]*4 #down*4/up*4
```

The uncertainty is asymmetric, so it uses the 8-element format of `ggH_higgsPt`. The first four values are the downward fractions for bins A, B, C, D and the last four are the upward fractions. The downward shift is 6.7% and the upward shift is 4.6%. The same value goes in every bin because it is a global cross section effect. `make_datacard` writes `1 - down` and `1 + up`, which gives `0.933/1.046` in the datacard.

Pitfalls: enter `0.067`, not `6.7`, and keep down before up.

### Blanks 2 to 4: observed counts in B, C, D

```python
b = np.count_nonzero(np.logical_and(Nhits[k] >= N_RECHIT_CUT, np.abs(deltaPhi[k]) >= DPHI_CUT))
c = np.count_nonzero(np.logical_and(Nhits[k] < N_RECHIT_CUT, np.abs(deltaPhi[k]) >= DPHI_CUT))
d = np.count_nonzero(np.logical_and(Nhits[k] < N_RECHIT_CUT, np.abs(deltaPhi[k]) < DPHI_CUT))
```

The bin layout is a 2x2 grid in (Nhits, deltaPhi):

| | dPhi < 0.75 | dPhi >= 0.75 |
|---|---|---|
| Nhits >= 130 | A (signal region) | B |
| Nhits < 130 | D | C |

The C and D assignment is forced by the ABCD line in `make_datacard`: `A rateParam chA bkg (@0*@2/@1) B,C,D`. In Combine, `@0`, `@1`, `@2` refer to B, C, D in the order listed, so A = B * D / C. That holds only if A/B = D/C, which requires D to be in the low-dPhi column and C in the high-dPhi column. The signal yields in the example datacard support this: A/B is about 2009/7.9 = 250 and D/C is about 286/0.84 = 340. The other pairing would give C/D around 0.003 for signal, which is inconsistent with A/B.

### Blanks 5 and 6: signal yields in C and D

```python
# bin C
cond = np.logical_and(Nhits[k]<N_RECHIT_CUT, deltaPhi[k]>=DPHI_CUT)
signal_rate.append(np.sum(w[cond]))

# bin D
cond = np.logical_and(Nhits[k]<N_RECHIT_CUT, deltaPhi[k]<DPHI_CUT)
signal_rate.append(np.sum(w[cond]))
```

The signal uses the same cuts as the data, so the signal and observation are in identical bins. `w` already includes the ctau reweighting, the cross section, and the luminosity, so the sum is the expected number of signal events in the bin.

Pitfalls: a mismatch between `>` and `>=` between data and signal, and swapping C and D (the datacard is still valid but the ABCD prediction for A is wrong).

Notebook markdown cells explaining these solutions were inserted after the corresponding code cells.

## Discussion 6.1: where is the ABCD relation constrained?

The B, C, D bins each have a free `rateParam` whose starting value is the observed count, with range `[0, 7x]`. The Poisson likelihood of the observed counts in B, C, D constrains those free yields. The bin A yield is not free. It is set by the formula `B * D / C`, so the observation in A is compared to this prediction. The ABCD relation is therefore constrained by the data in B, C, D, and the test in A is the signal region comparison.

## Question 6.1: what does `norm` do, and why scale the signal?

> Do you know what does the rateParam `norm` do and why do we add the rateParam to change the normalization of the signal yield?

The line `norm rateParam * signal 1` multiplies the signal rate in every bin by a free parameter `norm`, whose starting value is 1. The expected signal in bin i becomes `norm * sig_i`.

Reasoning:
1. The signal yield at BR = 1 varies from O(1) to O(1000) events in bin A across the ctau points, because the reconstruction efficiency depends strongly on the LLP decay position.
2. Combine's asymptotic limit is accurate when the fitted signal strength `r` is of order 1 to 15. If the expected signal is O(1000) events, the limit on `r` is around 0.001 and the approximation is poor. If the expected signal is 0.01 events, the limit on `r` is very large.
3. Setting `norm = 1/sigA` rescales the signal so that the expected signal in bin A is about 1 event. The limit on `r` then falls in the good range for every ctau.
4. Because the rescaling is a known constant, the physical limit is recovered by dividing the `r` limit by `sigA` (see the limit conversion below).

Why freeze it: if `norm` floats, the fit can absorb signal into the normalization and the limit on `r` no longer refers to a fixed signal model. Freezing it means `r` is the only signal parameter.

Note on the episode text: Solution 6.1 mentions `--setParameters norm=0.001`. That value is used for the background-only fit (it effectively removes the signal). For the limit calculation the value is `norm = 1/sigA`, as in `run_combine.py`.

## Results: background-only fit

```
combine -M MultiDimFit datacard.txt --saveWorkspace -n Snapshot --freezeParameters r,norm --setParameters r=0,norm=0.001
combine -M FitDiagnostics --snapshotName MultiDimFit --bypassFrequentistFit higgsCombineSnapshot.MultiDimFit.mH120.root --saveNormalizations --saveShapes --saveWithUncertainties
```

The first command makes a snapshot with `r = 0` and the signal effectively off. The second fits the background only, and `shapes_fit_b/chA/total_background` gives the predicted background in bin A. Compare its value and uncertainty with the observation in A. This is an interpretation step, not a code change.

## Run Higgs Combine (`scripts/run_combine.py`)

The blanks are the datacard directory and the freeze option.

```python
card_dir = CMSSW_BASE + "/src/MDS_CMSDAS/combine/datacards/"
```

The notebook writes datacards to `<repo>/combine/datacards/`, and the repo is at `$CMSSW_BASE/src/MDS_CMSDAS`. Because `output_tree_dir` is derived from `card_dir` by replacing `datacards` with `limitTrees`, the limit trees go to `<repo>/combine/limitTrees/`.

```python
comb_limit = "combine -M AsymptoticLimits  {card_dir}/{card}.txt -n.{card} --setParameters norm={norm} --freezeParameters norm"
```

`--freezeParameter=norm` was changed to `--freezeParameters norm`, which is the documented option name and matches the episode command. `norm` is set to `1/sigA`, where `sigA` is the first `rate` entry of the datacard (the signal yield in bin A at BR = 1).

Output: `higgsCombine.<name>.AsymptoticLimits.mH120.root`, moved into `combine/limitTrees/`.

## Make limit plots (`scripts/limitPlot.ipynb`)

Blanks:
- `dataCardDir = os.getcwd().replace('scripts', 'combine/datacards')`
- `limitDir = os.getcwd().replace('scripts', 'combine/limitTrees')`

These match the output locations above. The `h_obs` TGraph was already present in the skeleton, built from `limits[:,5]`, and matches the order of the six numbers (the observed limit is the last entry).

Limit conversion, the line marked `CHANGE THIS LINE`:
```python
limits[i] = np.array(T['limit'].array())/signal_norm
```
This is already correct and was kept. `signal_norm` is sigA at BR = 1, read from the datacard. Combine's `r` is the expected number of signal events in bin A (because `norm * sigA = 1` event), so `BR_lim = r_lim / sigA`.

Pitfall: multiplying instead of dividing gives an error of `sigA^2`.

## Discussion 6.1: shape of the limit

Short lifetimes: LLPs decay promptly, close to the interaction point, before the CSC. Few clusters are reconstructed in the CSC, so the signal efficiency is low and the limit on BR is weak.

Long lifetimes: LLPs decay beyond the detector or pass through it without decaying in the CSC volume. The efficiency again drops, so the limit weakens.

Between these, the decay probability in the CSC is largest and the limit is best.

## Exercise: 15 and 55 GeV

The notebooks and `run_combine.py` are written for the 40 GeV sample only. To repeat for another mass, change the following (assuming the same EOS naming pattern, which should be checked for the 15 and 55 GeV samples):
1. `create_datacard.ipynb`: `mass = [40]` to `[15]` or `[55]`, the `MS-40` sample path to the matching mass, and `mx40` in `modelName` if it is hard-coded.
2. `run_combine.py`: `sample_name = 'ggH_HToSSTobbbb_mh125_mx40'` to the new mass.
3. `limitPlot.ipynb`: `name = 'ggH_HToSSTobbbb_mh125_mx40'` to the new mass.
4. Plot the three mass curves on one canvas by adding a second and third set of TGraphs in `limitPlot.ipynb` with a different line colour and legend entry.

## Physics sanity checks

- Datacard format: `imax 4` (four counting bins). `jmax 1` (two processes, signal id 0 and background id 1). `kmax *` (number of nuisance lines is automatic). The `rate` line gives signal yields and 1 for the background (the background is set by the rateParams).
- Each nuisance line has two entries per bin (signal, background). A `-` means the nuisance does not apply to that process.
- Asymmetric lnN: 8 numbers (4 down, 4 up). Symmetric lnN: 4 numbers.
- `setParameters` sets the starting value of a parameter. `freezeParameters` removes it from the minimization, so it stays at its set value.
- Units: the signal yields are event counts for L = 137 fb^-1 and sigma = 48 pb at BR = 1. The ctau in filenames is in mm, and `limitPlot.ipynb` converts to m for the x axis. The y axis is BR(h -> SS), which is dimensionless.
