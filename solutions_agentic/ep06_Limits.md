# Episode 06 — Results and statistical analysis (datacards, Higgs Combine, limit plot)

This file documents the in-place solutions in `scripts/create_datacard.ipynb`, `scripts/run_combine.py` and `scripts/limitPlot.ipynb`; the datacard text itself is written by `helper.make_datacard()` in `lib/helper.py`. No code was executed to produce these answers - everything is derived from the repo (the example datacard printed in the notebook, the helper, and the variable/cut definitions of the earlier episodes).

Conventions used below (they come from the notebook and the helper, not from me):

| name | meaning |
| --- | --- |
| `Nhits` | `cscRechitClusterSize` of the leading selected cluster of the event |
| `deltaPhi` | `cscRechitClusterMet_dPhi`, angle between the cluster and the MET direction |
| `N_RECHIT_CUT = 130`, `DPHI_CUT = 0.75` | the ABCD boundaries (episode 05/06) |
| bin A | `Nhits >= 130` and `|deltaPhi| < 0.75` -> signal region |
| bin B | `Nhits >= 130` and `|deltaPhi| >= 0.75` -> sideband sharing `Nhits` |
| bin C | `Nhits < 130` and `|deltaPhi| >= 0.75` -> **crisscross** bin |
| bin D | `Nhits < 130` and `|deltaPhi| < 0.75` -> sideband sharing `deltaPhi` |
| `sigA` | signal yield in bin A, i.e. the second token of the `rate` line of a card |

With this labelling the ABCD product rule is `A = B x D / C`, which is exactly the formula the helper writes into the card (see section 3). Note one inconsistency in the exercise material itself: the "# ABCD validation" markdown of `scripts/ABCD_validation.ipynb` (episode 05) tells the student to compute `A = B x C / D`, i.e. it calls the crisscross bin D, while `lib/helper.py` calls it C. The two are the same physics with the labels C and D exchanged, but the datacards must use the definition above, because `helper.make_datacard()` hard-codes `(@0*@2/@1)` with the parameter list `B,C,D`. The ep05 solution fills the validation cell as `pred = b*d/c`, so the two notebooks agree with each other even though the ep05 prose does not.

---

## 1. Ex: add the QCD scale uncertainty (`create_datacard.ipynb`)

> **Can you fill in the QCD scale uncertainty using the same format as `ggH_higgsPt`?**

```python
for m in mass:
    for ct in OLD_CTAU:
            key = 'MC_'+str(m)+'_'+str(ct)

            ggH_higgsPt[key] = [0.20] * 4 + [0.13]*4 #down*4/up*4
            ggH_qcdScale[key] = [0.067] * 4 + [0.046]*4 #down*4/up*4 (6.7% down / 4.6% up, yellow report 4)
            ggH_pdf[key] = [0.032]*4
```

**Explanation.** The notebook markdown states the size of the uncertainty: the renormalisation/factorisation scale variation changes the ggH *total* cross section by 6.7% down and 4.6% up (LHC HXSWG Yellow Report 4). The asymmetric format used by `helper.make_datacard()` is a list of **8** numbers: `down` for bins A, B, C, D followed by `up` for bins A, B, C, D, and it prints `lnN` factors `1-down` / `1+up`. Mirroring `ggH_higgsPt` therefore gives `[0.067]*4 + [0.046]*4`, which produces `ggH_qcdScale lnN 0.933/1.046 -` - identical to the `ggH_qcdScale` line of the example card in the notebook, which is the check that the numbers and the order are right.

Pitfalls: (1) the entries are *relative* shifts and the downward one is entered as a positive fraction; (2) the list length must be 4 (symmetric) or 8 (asymmetric) or the helper hits `assert(False)`; (3) these uncertainties apply to the signal only - the `-` in the background column is what makes the background rate free instead, so it is constrained by the ABCD `rateParam`s rather than by a theory uncertainty.

---

## 2. Ex: complete the ingredients for the datacard (`create_datacard.ipynb`)

> Fill in everywhere that says `FILL IN ME` to complete the definition of bin B, C, D for data and then signal

Data (the unblinding), mirroring the given bin-A line and only flipping which side of each boundary is taken:

```python
k = 'data'
a = np.count_nonzero(np.logical_and(Nhits[k] >= N_RECHIT_CUT, np.abs(deltaPhi[k]) < DPHI_CUT))
b = np.count_nonzero(np.logical_and(Nhits[k] >= N_RECHIT_CUT, np.abs(deltaPhi[k]) >= DPHI_CUT))
c = np.count_nonzero(np.logical_and(Nhits[k] < N_RECHIT_CUT, np.abs(deltaPhi[k]) >= DPHI_CUT)) # C = crisscross
d = np.count_nonzero(np.logical_and(Nhits[k] < N_RECHIT_CUT, np.abs(deltaPhi[k]) < DPHI_CUT)) # D = low Nhits, small dphi
observation= np.array([a,b,c,d])
```

Signal yields inside the `(mass, ctau)` loop (bins A and B were given):

```python
 # bin A
cond = np.logical_and(Nhits[k]>=N_RECHIT_CUT, deltaPhi[k]<DPHI_CUT)
signal_rate.append(np.sum(w[cond]))

# bin B
cond = np.logical_and(Nhits[k]>=N_RECHIT_CUT, deltaPhi[k]>=DPHI_CUT)
signal_rate.append(np.sum(w[cond]))

# bin C : the crisscross corner, low Nhits and LARGE dphi
cond = np.logical_and(Nhits[k]<N_RECHIT_CUT, deltaPhi[k]>=DPHI_CUT)
signal_rate.append(np.sum(w[cond]))

# bin D : low Nhits and SMALL dphi
cond = np.logical_and(Nhits[k]<N_RECHIT_CUT, deltaPhi[k]<DPHI_CUT)
signal_rate.append(np.sum(w[cond]))
```

**Explanation.** Both lists must have exactly four entries because `make_datacard()` asserts `len(signal_rate) == len(observation) == 4` and writes one column per channel (`chA`...`chD`). For data the yield is a simple count (`np.count_nonzero`); for MC it is the **sum of the event weights** `w = weight * weight_ctau`, because `weight = L x sigma / N_gen` (with the pileup and MET scale factors folded in) and `weight_ctau = helper.weight_calc(...)` reweights the sample from its generated lifetime to the requested one. Counting MC events instead of summing weights would return the yield for the generated luminosity/lifetime.

How do we know which corner is C and which is D? The example card printed in the notebook is self-consistent only with C being the crisscross bin:

* the observation line is `3, 3, 121, 50`, so bin C = 121 is by far the most populated corner - that is the anti-signal region (small clusters, back-to-back with the MET), and `50 : 121 ~ 1 : 2.4` is also about the `0.75 : (pi - 0.75)` phase-space ratio expected for a background that is uniform in `|dphi|`;
* the signal rates `2.0098e+03, 7.889, 0.8444, 286.0` satisfy the independence of the two variables: `sigC ~ sigB x sigD / sigA = 7.889 x 286.0 / 2009.8 = 1.1` versus the actual `0.84`, while the alternative labelling would need `sigD ~ 0.003` instead of `286`.

Pitfalls: (1) the signal masks use `deltaPhi[k] < DPHI_CUT` without `np.abs`, exactly as the given bin-A/bin-B lines do, while the data lines use `np.abs(...)` - keep each style where it already was used rather than "fixing" one of them; (2) `observation` is computed **once** outside the loop (background does not depend on the signal point), while `signal_rate` is rebuilt for every `ctau`; (3) a control bin with zero counts makes the flat prior `[0, 0]` and Combine will fail.

---

## 3. Discussion 6.1 - "Do you understand the example datacards? Where is the ABCD relationship constrained for background?"

```
imax 4
jmax 1
kmax *
shapes * * FAKE

bin      chA     chB     chC     chD
observation        3.00            3.00          121.00           50.00

bin      chA     chA     chB     chB     chC     chC     chD     chD
process          signal          bkg     signal          bkg     signal          bkg     signal          bkg
process          0       1       0       1       0       1       0       1
rate     2.009794e+03    1       7.888562e+00    1       8.444456e-01    1       2.860445e+02    1

A        rateParam       chA     bkg      (@0*@2/@1)                    B,C,D
B        rateParam       chB     bkg     3.00        [0,21.00]
C        rateParam       chC     bkg     121.00        [0,847.00]
D        rateParam       chD     bkg     50.00        [0,350.00]

ggH_higgsPt      lnN     0.792/1.134     -       ...
ggH_qcdScale     lnN     0.933/1.046     -       ...
ggH_pdf          lnN     1.032   -       ...

norm rateParam * signal 1
```

**Explanation, line by line.**

* `imax 4 / jmax 1 / kmax *`: four bins (`chA`...`chD`), one background process (`jmax` counts the non-signal processes), an unlimited number of nuisance parameters. `shapes * * FAKE` says there are no templates, so this is a four-bin counting experiment.
* `observation`: the unblinded data counts in the four corners, the numbers `a, b, c, d` of section 2.
* `bin`/`process`/`process`: eight columns, a signal and a background one per bin. The second `process` row (`0 1 0 1 ...`) is what tells Combine that the signal-strength modifier `r` multiplies only the `0` (signal) columns.
* `rate`: the expected yields. The signal entries are the weighted MC yields of section 2 (already normalised to `L x sigma / N_gen`, i.e. to BR(h -> SS) = 1). The background entries are the dummy `1`: they are overwritten by the `rateParam` lines below, so the background column carries no external prediction at all.
* `ggH_*  lnN`: signal-rate theory uncertainties (section 1), `-` in every background column.

**The ABCD constraint lives in these four lines:**

```
A   rateParam   chA   bkg   (@0*@2/@1)   B,C,D
B   rateParam   chB   bkg   3.00         [0,21.00]
C   rateParam   chC   bkg   121.00       [0,847.00]
D   rateParam   chD   bkg   50.00        [0,350.00]
```

`B`, `C` and `D` are nuisance parameters: free normalisations of the `bkg` process in their own channel, with nominal values taken from the observed counts of the three control regions and flat priors `[0, 7 x observation]` (so nothing but the data pulls on them, and the upper bound keeps the fit from running away). The first line is the ABCD method: the background yield in the signal region is *not* a free parameter, it is the expression `(@0*@2/@1)` evaluated on the parameters named at the end of the line, i.e. `@0=B`, `@1=C`, `@2=D`, giving

```
n_bkg(chA) = B x D / C        (example: 3 x 50 / 121 = 1.24)
```

So the fit "measures" the three control regions from data and *predicts* the signal region from the product rule, inside a single binned Poisson likelihood. That is the advantage over the closure test of episode 05: the uncertainty on the prediction comes from the Poisson terms of all four bins automatically (its Gaussian limit is `(B x D / C) x sqrt(1/B + 1/C + 1/D)`, the `unc_pred` the students wrote in the ABCD validation notebook).

Pitfalls: the `@` indices follow the **order of the names on the line**, so reordering `B,C,D` silently changes `B*D/C` into another combination; and a control region with zero observed events produces the degenerate prior `[0, 0]`.

---

## 4. Question 6.1 - "What does the rateParam `norm` do and why do we add it?"

```
norm rateParam * signal 1
```

```bash
combine -M AsymptoticLimits card.txt -n.TAG --freezeParameters norm --setParameters norm=0.000497
# 0.000497 = 1./sigA for the example card (sigA = 2.009794e+03)
```

**Explanation.** `norm` multiplies the signal rate in *all* channels (the `*` in the third column is the wildcard for the bins), so the number of signal events used by the fit is `r x norm x sig_i`, where `r` is the signal strength Combine reports. In principle `norm` could float - then it would simply absorb any normalisation uncertainty on the signal - but here it is deliberately **frozen** (`--freezeParameters norm`) and set to a fixed value (`--setParameters norm=...`).

The reason is numerical/statistical, not physical: the signal yield of this analysis varies from O(1) to O(1000) events across the scanned LLP lifetimes (see the `rate` lines of the cards: `sigA` alone spans that range), while the signed profile-likelihood-ratio machinery behind `-M AsymptoticLimits` is an asymptotic approximation that works best when the fitted signal strength is of order 1-15. The trick is therefore to scale the signal template to O(1) by setting `norm = 1./sigA` (the "rule of thumb" quoted in the episode) and to let `r` absorb the remaining signal-strength factor. Since the MC yields already assume BR(h -> SS) = 1, the physical limit is

```
BR(h -> SS) = r x norm = r / sigA
```

which is exactly the multiplication applied in `limitPlot.ipynb` (`limits[i] = np.array(T['limit'].array())/signal_norm`). Pitfall: forgetting to multiply the scale factor back gives limits that are wrong by up to three orders of magnitude, and the size of the error depends on ctau, so the plot would still "look" like a limit curve.

---

## 5. Background-only fit: does the observation agree with the ABCD prediction?

```bash
combine -M MultiDimFit datacard.txt --saveWorkspace -n Snapshot --freezeParameters r,norm --setParameters r=0,norm=0.001
combine -M FitDiagnostics  --snapshotName MultiDimFit --bypassFrequentistFit higgsCombineSnapshot.MultiDimFit.mH120.root --saveNormalizations --saveShapes --saveWithUncertainties
```

```bash
root -l fitDiagnosticsTest.root
shapes_fit_b->cd()
chA->cd()
total_background->GetBinContent(1) # this will print the background prediction
total_background->GetBinError(1) # this will print the uncertainty on the background prediction
```

**Explanation.** The first command performs a `MultiDimFit` with the signal strength pinned to zero (`--freezeParameters r,norm --setParameters r=0,norm=0.001`), i.e. a pure background-only fit, and saves the resulting workspace with `-n Snapshot` (hence the file name `higgsCombineSnapshot.MultiDimFit.mH120.root`). The second command runs `FitDiagnostics` on that saved snapshot instead of refitting (`--snapshotName MultiDimFit --bypassFrequentistFit`) and asks it to dump the post-fit yields per bin with their uncertainties (`--saveNormalizations --saveShapes --saveWithUncertainties`), producing `fitDiagnosticsTest.root`.

Inside `shapes_fit_b/chA`, `total_background` is the post-fit background prediction in the signal region: `GetBinContent(1)` is the central value, `GetBinError(1)` the total (statistical + nuisance-propagated) uncertainty. Compare it with the `observation` value of bin A in the card: agreement within the error means the ABCD model describes the data and a signal region excess, if any, does not have to be explained by a bad background estimate. The fitted `B`, `C`, `D` normalisations in `shapes_fit_b`/the `normalizations` tree should also sit close to their observed counts, since the pull on them is only Poisson.

Pitfall: `FitDiagnostics` writes `fitDiagnosticsTest*.root` and `MultiDimFit` writes `higgsCombine*` in the *current* directory, so run both commands in the same folder (or move the snapshot file) or the second command cannot find its input.

---

## 6. `scripts/run_combine.py` - the two filled spots

The datacard directory (`create_datacard.ipynb` writes the cards with `outDataCardsDir = os.getcwd().replace('scripts', 'combine/datacards/')`):

```python
    card_dir= CMSSW_BASE + "/src/MDS_CMSDAS/combine/datacards/"

    output_tree_dir = card_dir.replace('datacards', 'limitTrees')
```

and the Combine invocation:

```python
        comb_limit = "combine -M AsymptoticLimits {card_dir}{card}.txt -n.{card} --setParameters norm={norm} --freezeParameters norm".format(card_dir=card_dir,card=name,norm = 1./norm)
```

**Explanation.** `os.getenv('CMSSW_BASE')` has no trailing slash, so the appended string must start with `/`; the trailing `/` is kept because the existence check and the `readNorm()` call concatenate `card_dir + name + '.txt'` directly, and the directory name must literally contain `datacards` for the `.replace('datacards', 'limitTrees')` trick that creates `$CMSSW_BASE/src/MDS_CMSDAS/combine/limitTrees/`. The card name built in the loop, `ggH_HToSSTobbbb_mh125_mx40_ctau<ct>mm`, matches the `modelName` of the notebook exactly, so every lifetime point is picked up; `readNorm()` returns `sigA` from the card and the script runs Combine with `norm = 1./sigA`, the scaling of section 4. With the tag `-n.{card}` Combine writes `higgsCombine.<card>.AsymptoticLimits.mH120.root`, which is what the final `mv` moves into `limitTrees/` (and what `limitPlot.ipynb` opens).

Two things were fixed in the original line: `--freezeParameter=norm` is not a Combine option (the option is the plural `--freezeParameters`, with the names separated by spaces), and `{card_dir}/{card}.txt` would have produced a doubled slash once `card_dir` ends in `/`.

Pitfalls: `--setParameters norm={norm}` prints the float in `repr` form (e.g. `0.0004975625...`), which is what Combine expects, but it also means the value used for the plot must come from the same card - that is why `limitPlot.ipynb` re-reads `sigA` with its own copy of `readNorm()` instead of hard-coding it. The script's `--dryRun` option is parsed but never used; uncomment `#print(comb_limit)` if you want to inspect the commands first.

---

## 7. `limitPlot.ipynb` - the limit trees, the observed graph and the r -> BR conversion

Filled in (the two placeholder paths of the loading cell):

```python
dataCardDir = os.getenv('CMSSW_BASE') + '/src/MDS_CMSDAS/combine/datacards'
limitDir = os.getenv('CMSSW_BASE') + '/src/MDS_CMSDAS/combine/limitTrees'
```

The rest of the pipeline, with the line the notebook marks `CHANGE THIS LINE` and the `h_obs` graph the markdown asks for:

```python
        T = uproot.open(limitTrees[m])['limit']
        limits[i] = np.array(T['limit'].array())/signal_norm #CHANGE THIS LINE
...
h = rt.TGraph(len(ctaus_m[cond]),ctaus_m[cond], limits[:,2][cond])
h_obs = rt.TGraph(len(ctaus_m[cond]),ctaus_m[cond], limits[:,5][cond])
h_exp1sig = rt.TGraph(2*len(ctaus_m[cond]), np.hstack((ctaus_m[cond], np.flip(ctaus_m[cond]))), np.hstack((limits[:,1][cond], np.flip(limits[:,3][cond]))))
```

**Explanation.** Each `higgsCombine.<card>.AsymptoticLimits.mH120.root` contains a `limit` tree with six rows in the order `-2 sigma, -1 sigma, median expected, +1 sigma, +2 sigma, observed` - which is why `limits` has shape `(n_ctau, 6)`. `h_obs` therefore takes **column 5** (the last row, the observed limit) exactly where `h` (the median expected) takes column 2, and the +-1 sigma band takes columns 1 and 3 with the ctau axis walked forwards and then backwards so the `TGraph` closes. Note that the expected/observed distinction is *positional*, not a selection on `q0`: it comes from the row order in which `AsymptoticLimits` fills the tree. The conversion to a branching ratio is the division by `signal_norm = sigA` (= multiplication by the frozen `norm`), so nothing more has to be done for `h_obs`: it is on the same BR scale as `h`.

In the notebook as shipped, `h_obs` was already written and the two paths were the only blanks; the exercise code for the 15/55 GeV overlay was added as the last cell.

Pitfalls: unfilled points keep the sentinel `100000` and are dropped by `cond = limits[:,2] < 1000000`, so a curve that stops early means Combine failed (or `sigA` was 0) for those lifetimes; and the loading cell checks `os.path.isfile` only for the *card*, so a missing limit tree raises - add `os.path.isfile(limitTrees[m])` if you expect incomplete runs. Also note the ctau list in this notebook (10 mm to 10 m, 16 points) is a subset of the 40 points that `run_combine.py` runs, so the extra cards are simply never read.

---

## 8. Discussion 6.1 (plot) - "Did you expect the shape of the limit to look like this? Why does the sensitivity decrease for very long or very short lifetimes?"

**Explanation.** The curve is an inverted umbrella: weak limits at both ends of the ctau axis and a best sensitivity in the middle (for this analysis, of order 10 cm to a few m), and the same shape is obtained for every mS.

The reason is the probability that the LLP decays *inside the instrumented muon system and in front of the clusters we call signal*. In the lab frame the decay length is `L = c*tau x beta*gamma`, and for a decay between radius `R_min` and `R_max` the probability is

```
P(R_min < L < R_max) = exp(-R_min/(beta*gamma*c*tau)) - exp(-R_max/(beta*gamma*c*tau))
```

which rises from zero, peaks when `beta*gamma*c*tau` is comparable to the size of the muon system, and falls again as `~ 1/c*tau` once the LLP is stable on detector scales. At very short ctau (mm and below) the LLP decays at the interaction point or in the calorimeters: the b-jets look ordinary, nothing reaches the CSCs, and the signal efficiency of the MDS cluster selection collapses. At very long ctau (km) the LLP escapes CMS without decaying, again leaving no cluster. In between, essentially every h -> SS decay that points into the acceptance produces a large, delayed, MET-aligned cluster, which is why the limit is strongest there. The y axis is a limit on BR(h -> SS), so the shape of the curve is a pure acceptance/efficiency statement multiplied by 1/(expected signal yield).

At fixed ctau the curve also moves with mS: `p* = sqrt(mH^2/4 - mS^2)` is the momentum of each S in the Higgs rest frame, so a heavier S is slower and its `beta*gamma` (hence its lab decay length and where it decays) is smaller, while a very light S is strongly boosted and its two b-jets and their hits are collimated, which changes the cluster size and the efficiency of the cut-based ID. The student should check the plot against slide 25 of the introduction slides by reading the y value at the minimum of their own curve, and check the axis convention first (`c#tau [m]`, log-log) - the repo does not state a number to compare to.

---

## 9. Exercise - "create datacards, run combine and make limits plot for 15 and 55 GeV ... plot them in the same canvas"

```bash
# 1) create_datacard.ipynb: extend the mass list
mass = [15, 40, 55]
#    and make the input file name follow m (it is hard coded to MS-40 today):
#    fpath[key] = path + 'ggH_HToSSTobbbb_MH-125_MS-'+str(m)+'_ctau-'+str(ct)+'_137000pb_weighted.root'
#    then re-run all cells (the observation does not change, only signal_rate does)

# 2) run_combine.py, once per mass (the cards and trees of all masses share the same two folders)
python run_combine.py     # with sample_name = 'ggH_HToSSTobbbb_mh125_mx15'
python run_combine.py     # with sample_name = 'ggH_HToSSTobbbb_mh125_mx40'
python run_combine.py     # with sample_name = 'ggH_HToSSTobbbb_mh125_mx55'
```

and, in `limitPlot.ipynb`, the last cell added to the notebook builds one observed `TGraph` per mass from `masses = [15, 40, 55]` and draws them on a single `TCanvas` (first with `LPA` to define the frame, the others with `Lsame`).

**Explanation.** Three things have to be generalised, and they are all in the file names, not in the physics: the mass list of the loading cell of `create_datacard.ipynb` (plus the hard-coded `MS-40` substring of the input file name, otherwise the same mS = 40 ntuples are used for all three "masses" and the cards would differ only in name), the `sample_name` of `run_combine.py`, and the card/tree name prefix in the plot notebook. The ABCD observation, the theory uncertainties and the ctau list stay identical, so the three curves differ only through the signal yields - which is exactly the comparison the exercise is about.

Pitfalls: check that the `MS-15` / `MS-55` weighted ntuples actually exist under the EOS path before assuming a blank plot; and verify that the generated lifetimes `OLD_CTAU` cover the reweighting range you ask for (`weight_calc` extrapolates the shape with an exponential, which is only trustworthy inside the generated window). If you re-run `create_datacard.ipynb` with a longer mass list, remember that it overwrites cards in place - delete `combine/datacards/` first if you want to be sure that no stale card from an earlier version is picked up by `run_combine.py`.

---

## 10. Notes for the instructor (found, not changed)

* **Time-window mismatch (needs instructor decision):** `create_datacard.ipynb` requires `-5 ns < cscRechitClusterTime < 12.5 ns` for its ABCD/SR counts, but the episode-04 text defines the signal region as "cluster time between -5 ns and -12.5 ns" (an *early* window), and `ABCD_validation.ipynb` uses `t < -12.5` ns for the OOT validation region, consistent with the episode text. The datacard window `(-5, +12.5)` includes the in-time peak and is not adjacent to the OOT region, so the three definitions cannot all be right; the sign in the datacard notebook (and in one heading of `analysis_strategy.ipynb`) looks like a material typo. Not changed (pre-existing skeleton code, outside the blanks).
* The example card printed in `create_datacard.ipynb` shows `ggH_higgsPt lnN 0.792/1.134`, while the code cell hard-codes 0.20 down / 0.13 up (i.e. `0.80/1.13`); the printed example is illustrative. The `ggH_qcdScale 0.933/1.046` line of the example is the one that pins down the convention (1-0.067 / 1+0.046).
* `helper.make_datacard()` writes the four corner bins as `chA...chD` with the constraint `(@0*@2/@1)` = B x D / C, i.e. it treats **C** as the crisscross bin; the prose of `scripts/ABCD_validation.ipynb` says `A = B x C / D`, which treats **D** as the crisscross bin. The notebook skeletons disagree; the filled solutions of episodes 05 and 06 both use C = crisscross so that the cards are consistent with the helper. Worth an instructor decision.
* `run_combine.py` defines a `--dryRun` option that is never used, and an `insert()` helper that is never called.
* `event_display.ipynb` still reads the `/2024/long_exercises/` EOS path while `create_datacard.ipynb` uses `/2025/long_exercises/`.
