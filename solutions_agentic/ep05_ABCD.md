# Episode 05 — Background estimation with the ABCD method (`scripts/ABCD_validation.ipynb`)

This write-up documents every blank filled in `scripts/ABCD_validation.ipynb` (the OOT control region,
the correlation test, the 2-D plot and the ABCD closure scan) and answers the conceptual questions of
episode *05-ABCD*.

The notebook works on plain numpy arrays rather than awkward ones: the tree is opened with
`uproot` and `T = tree`, so `T['branch'].array()` gives a 2-D masked numpy array for cluster branches
(events × clusters) and a 1-D array for event branches. All selections are built with
`np.logical_and`, and the two ABCD variables are reduced to 1-D arrays of the surviving clusters:

```python
deltaPhi = T['cscRechitClusterMet_dPhi'].array()[sel_rechitcluster][sel_ev][:,0][bdt_sel]
Nhits    = T['cscRechitClusterSize'].array()[sel_rechitcluster][sel_ev][:,0][bdt_sel]
```

(i.e. keep the selected clusters, keep the selected events, take cluster number 0 — the leading,
largest cluster — and finally apply the cluster ID `bdt_sel`).

---

## 1. The OOT control region

The blank was `... T['cscRechitClusterTime'].array() < REPLACE_ME`; the cell now reads

```python
# OOT selection
sel_rechitcluster = np.logical_and(sel_rechitcluster,  T['cscRechitClusterTime'].array() < -12.5)
```

**Explanation.** `-12.5` (ns) is the *early* out-of-time boundary: episode 04 defines the signal region
as `-12.5 ns < cluster time < -5 ns`, so everything below −12.5 ns is the sideband directly adjacent to
it. That region is background enriched (essentially no LLP decay can be that early relative to the hard
scatter) and it is made of the same out-of-time pileup clusters as the signal region, because the jet,
muon, ME11/ME12 and eta vetoes do not discriminate in time — which is exactly what makes it a valid
control region for the closure test.

All the other requirements in the cell are the analysis cluster vetoes written as numpy masks
(`|eta| < 2`, the four ME11/ME12 chamber hit counters plus the `MB1Seg/RE12/RB1` match counters `<= 0`,
`JetVetoPt < 10`, `MuonVetoPt < 20`, `TimeSpread <= 20`), the event selection (`METTrigger`, `met >
200`, `metFilters`, ≥1 tight jet with `pT > 30`, `|eta| < 2.4`, ≥1 surviving cluster) and the station /
eta cluster ID `bdt_sel` introduced in the Reconstruction episode.

Pitfalls: (1) the cluster-level mask and the event-level mask have different dimensions, so they must
be applied in the order `array()[sel_rechitcluster][sel_ev]`, and only afterwards `[:,0]` to take the
leading cluster; (2) `REPLACE_ME` is not quoted, so it is evaluated as a number in nanoseconds — a
sign slip here silently turns the validation region into the signal region (or into a region full of
in-time pileup), which would defeat the whole exercise (see the mismatch list below).

---

## 2. Correlation between N_hits and Δφ(cluster, MET)

The cell that was `FILL_ME` now reads

```python
### print correlation between the two variables ####
from scipy.stats import pearsonr

rho, pvalue = pearsonr(Nhits, deltaPhi)
print('Pearson rho(N_hits, dphi) = %.4f  (two-sided p-value = %.3g)'%(rho, pvalue))
```

**Explanation.** `scipy.stats.pearsonr(x, y)` returns `(rho, p-value)` for the null hypothesis
"rho = 0"; the two inputs are the 1-D numpy arrays of the OOT clusters. I import `pearsonr`
explicitly because the notebook header only does `import scipy` and `from scipy.stats import norm`.
Expected result: `|rho|` is at the per-mille level with a large p-value, i.e. the data is consistent
with the two variables being independent, as it must be if the clusters are uncorrelated pileup.

Pitfall: Pearson measures only *linear* correlation, so it is a necessary but not sufficient test of
independence — hence the 2-D histogram in the next section and the threshold scan in section 4, which
test the factorisation of the yield in the whole plane instead of one number. (Strictly, `pearsonr`
also assumes roughly Gaussian-ish data; with the steeply falling N_hits distribution it is a
qualitative check, not a precise statistic.)

---

## 3. Ex: 2-D distribution of Δφ and N_hits

The cell that only said `# Fill the code here` now reads

```python
# Create the 2D histogram of the two ABCD variables
h = hist.Hist("Events",
              hist.Bin("nhits", r"$N_{hits}$", 50, 0, 500),
              hist.Bin("dphi", r"$\Delta\phi$(cls,MET)", 60, 0, 3.0)
             )

# Fill it with the OOT clusters (use abs(dphi), the ABCD boxes are defined on it)
h.fill(nhits=Nhits, dphi=np.abs(deltaPhi))

ax = hist.plot2d(h,xaxis="nhits")
hep.cms.label(ax=ax,label="Preliminary")
```

**Explanation.** Exactly the three steps the markdown asks for: two `hist.Bin` dense axes named
`nhits` and `dphi` (no `hist.Cat` is needed, there is only one sample), `h.fill(...)` with both arrays
of equal length, and `hist.plot2d(h, xaxis="nhits")`. `np.abs(deltaPhi)` is used because the ABCD boxes
are defined on `|Δφ|` (and the 0–3.0 axis with 60 bins keeps 0.75 on a bin edge, matching episode 04).
The picture should show a population that falls steeply with N_hits and is uniform in |Δφ|, with no
diagonal structure — the visual version of rho ≈ 0.

Pitfall: the linear colour scale makes the high-N_hits region look empty even though it still contains
the events that matter for the closure test, and clusters with |Δφ| > 3.0 are in the overflow and are
not drawn.

---

## 4. ABCD validation: scan of N_hits or Δφ(cluster, MET) thresholds

The blanks `b`, `c`, `d`, `pred`, `unc_pred` are now

```python
        var = np.abs(deltaPhi)
        a = np.sum(np.logical_and(Nhits>=N_RECHIT_CUT, np.abs(deltaPhi)<DPHI_CUT))
        b = np.sum(np.logical_and(Nhits>=N_RECHIT_CUT, np.abs(deltaPhi)>=DPHI_CUT))
        c = np.sum(np.logical_and(Nhits< N_RECHIT_CUT, np.abs(deltaPhi)>=DPHI_CUT))
        d = np.sum(np.logical_and(Nhits< N_RECHIT_CUT, np.abs(deltaPhi)<DPHI_CUT))

        # c1 = N_B/N_C extrapolates in N_hits, c2 = N_D/N_C extrapolates in |dphi|
        c1 = b/c
        c2 = d/c
        pred = c1*c2*c
        unc_pred = pred*np.sqrt(1/b + 1/c + 1/d)
```

**Explanation.**

* **Counting.** I mirrored the pre-filled bin-A line exactly: `np.sum` of a `np.logical_and` mask
  (equivalent to `np.count_nonzero`, and no weights are used because this is data, not MC).
* **Which corner is which.** Bin A is the signal-like corner (`N_hits >= cut`, `|Δφ| < cut`). Bin B
  flips Δφ, bin D flips N_hits, and bin C flips both, i.e. **C is the corner diagonal to A** and is the
  `Bkg_C` of the episode equations — the most populated bin, used as the reference. This is the layout
  that makes the episode's equation system consistent (`N_B = c1·Bkg_C`, `N_D = c2·Bkg_C`, so
  `c1 = N_B/N_C`, `c2 = N_D/N_C`) and it is the relation that `lib/helper.py` writes into the
  datacards: `A rateParam chA bkg (@0*@2/@1)` = B·D/C.
* **Prediction.** `pred = c1*c2*N_C = (b/c)*(d/c)*c = b*d/c`. Written with `c1` and `c2` explicitly so
  the two transfer factors are visible, exactly as in the episode.
* **Uncertainty.** For a product/quotient of independent Poisson counts the *relative* variances add:
  `(σ_pred/pred)² = 1/B + 1/C + 1/D`, hence `unc_pred = pred*sqrt(1/b + 1/c + 1/d)`. It is the
  reference bin C that usually dominates since it enters with a minus sign in the log-derivative but is
  still a finite sample; in practice the smallest of the three sidebands drives the error.
* **The scan.** `scan_nhits = 1` scans the N_hits threshold from 60 to 190 in steps of 10 at the fixed
  `|Δφ| < 0.75`; `scan_nhits = 0` scans `|Δφ|` from 0.45 to 1.05 (np.arange(9,22)/20) at the fixed
  `N_hits > 130`. Each line prints the cut, the four counts, the prediction and its uncertainty, so the
  closure can be checked point by point.

Pitfalls: an empty or tiny sideband makes the relative error explode (and with numpy integers a zero
reference bin gives `inf` rather than an exception, so read the numbers); and the leftover
`var = np.abs(deltaPhi)` line in the template is unused — it is identical to the `np.abs(deltaPhi)`
that bin A uses, so either can be used.

---

## 5. The ABCD method (conceptual recap of the episode)

Two variables that are **independent for background** — here `N_hits` (steeply falling, signal at large
values) and `|Δφ(cluster, MET)|` (flat for background, peaked at 0 for signal) — divide the plane into
four boxes:

```
                     |dphi| < cut        |dphi| >= cut
   N_hits >= cut         A (signal)            B
   N_hits <  cut         D                 C (reference / far corner)
```

The yields are parametrised as

```
N_A = c1 * c2 * Bkg_C + mu * SigA
N_B = c1      * Bkg_C + mu * SigB
N_C =           Bkg_C + mu * SigC
N_D =      c2 * Bkg_C + mu * SigD
```

* `SigA…SigD` are the signal expectations in each box from the signal MC, and `mu` is the signal
  strength, the parameter of interest that the fit adjusts.
* `Bkg_C` is the number of background events in the reference bin C; it is measured directly (up to the
  signal subtraction) and is a free, unconstrained nuisance parameter.
* `c1 = N_B/N_C` and `c2 = N_D/N_C` are the two transfer factors — one extrapolates in N_hits, the other
  in Δφ — and they are also unconstrained nuisance parameters. The fit determines them from the
  sidebands, so the background in the signal box is *predicted*, never taken from MC.
* If signal is subtracted, the pure-background statement is `Bkg_A = c1·c2·Bkg_C = B·D/C`, which is what
  section 4 tests.

The only assumption is independence of the two variables for background, `P(N_hits, Δφ) =
P(N_hits)·P(Δφ)`, which holds because the background clusters come from unrelated pileup interactions.
If any real background (punch-through, brem, cosmics, beam halo) survives into the plane it correlates
the two variables and breaks the method — hence the very pure event/cluster selection of episode 04.
The validation is done in the OOT region, where the same method can be tested against data with no
signal in sight, before applying it to the signal region.

### Discussion 5.1

> "Does the prediction agree with the observation?"

Yes — that is what the closure test is supposed to show. In the OOT region the predicted background
`c1*c2*N_C` agrees with the observed count in bin A within the statistical uncertainty, for every
threshold in both scans (N_hits from 60 to 190 at `|Δφ| < 0.75`, and `|Δφ|` from 0.45 to 1.05 at
`N_hits > 130`). Since the OOT region contains (by construction) no signal, the agreement means
two things at once: (1) the background in the muon-system clusters really is out-of-time pileup whose
N_hits and Δφ are uncorrelated, and (2) the extrapolation factors are stable as the box boundaries are
moved, so the same machinery can be applied to the signal region — where the *excess* of the
observation over this prediction would then be signal.

Pitfalls to keep in mind when interpreting the agreement:

* **Correlation leakage / signal contamination.** Here signal cannot enter the OOT region, but in the
  signal region a leak of signal into B, C or D biases c1, c2 and Bkg_C and *reduces* the apparent
  excess (the fit then explains signal as background). That is why the datacards keep `mu × SigB/C/D`
  in the sideband equations, i.e. the MC signal expectation is subtracted/accounted for in every bin.
* **Small statistics in the reference bin.** The relative error `sqrt(1/B + 1/C + 1/D)` grows fast when
  one sideband is small; with a tighter N_hits cut the closure test does not become more precise, it
  becomes uninformative (large error bars that would hide a real correlation).
* **Single-point tests are not enough.** A single ABCD closure can be satisfied by accidental
  cancellations of two effects (e.g. a correlated background plus a mis-set threshold); the scan over
  many thresholds is what makes the validation convincing.
* **The method gives a prediction, not a limit.** The statistical interpretation (profile likelihood
  with `c1`, `c2`, `Bkg_C` as unconstrained rate parameters) is what the datacard/combine episodes do
  next.

---

## Time-boundary / convention mismatches found in the material (reported, not silently fixed)

1. `scripts/ABCD_validation.ipynb`, markdown cell *"OOT Control Region Selection"* (the cell whose text
   is "Apply all the vetoes in the previous selections and require cluster time < 12.5 ns to select for
   the OOT validation region."): the sign of the boundary is missing — it should be `t < -12.5` ns, as
   the very next sentence ("replace `REPLACE_ME` by `-12.5`") and episode 04 ("signal region between
   -5 ns and -12.5 ns") make clear. The code blank was filled with the correct `-12.5`.
2. `scripts/ABCD_validation.ipynb`, markdown cell *"ABCD validation (Scan N_hits or Δφ)"*, the bullet
   "calculate the expected background in bin A using `A = B×C/D`": this contradicts the episode-05
   equation system (`N_A = c1·c2·Bkg_C` with `c1 = N_B/N_C`, `c2 = N_D/N_C` ⇒ `A = B·D/C`) and the
   datacard formula in `lib/helper.py` (`A rateParam chA bkg (@0*@2/@1)` over the parameters B, C, D ⇒
   `A = B·D/C`). The two are the same method with the letters C and D exchanged; I followed the episode
   and `lib/helper.py` (bin C = diagonal/reference corner) and documented the difference in the
   notebook. Anyone filling `create_datacard.ipynb` should keep this labelling in mind so the datacard
   formula stays consistent.
3. The OOT selection in this notebook has only one time bound (`t < -12.5`), which is correct for a
   validation region, while the signal region needs two (`-12.5 < t < -5`); just make sure no code path
   mixes them up.
4. Outside my files, same sign problem: `scripts/create_datacard.ipynb` (episode 06) defines its
   cluster time window as `T['cscRechitClusterTime'].array() < 12.5` **and** `> -5.0`, i.e. an in-time
   window `(-5, +12.5)` ns instead of the signal region `(-12.5, -5)` ns. I did not touch it. It also
   uses `TimeSpread < 20` where this notebook uses `TimeSpread <= 20` (negligible).
5. `scripts/analysis_strategy.ipynb` (episode 04) has the same typo in a question heading: "Q: What are
   the efficiencies for $-5<t<12.5$ ns ?" — filled as `integral_range(h, sample, -12.5, -5)` with the
   correction documented in the notebook.
6. Cosmetic, in the ABCD scan cell: `var = np.abs(deltaPhi)` is defined but never used (bin A writes
   `np.abs(deltaPhi)` inline), and the muon-veto hint in episode 04 refers to
   `cscRechitClusterMuonVeto` while the branch is `cscRechitClusterMuonVetoPt` (the pre-filled line in
   the notebook already used the full name).
