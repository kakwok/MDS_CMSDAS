# Episode 05: ABCD background estimate (agentic solution)

Notebook: `scripts/ABCD_validation.ipynb`. The question file `_episodes/05-ABCD.md` is not in this repository or on the `2025` branch, so the questions below are quoted from the notebook markdown cells. The code was not executed, so no numbers are quoted.

Units: cluster time in ns, veto pT in GeV, |dphi| in radians, N_hits dimensionless.

---

## OOT control region selection (REPLACE_ME)

Question: "Apply all the vetoes ... and require cluster time < 12.5 ns to select for the OOT validation region. You will replace the `REPLACE_ME` variable by `-12.5`."

```python
sel_rechitcluster = np.logical_and(sel_rechitcluster,  T['cscRechitClusterTime'].array() < -12.5)
```

Explanation: the time cut keeps clusters below -12.5 ns, the out-of-time side of the in-time window. The vetoes (jet veto < 10 GeV, muon veto < 20 GeV, |eta| < 2, time spread <= 20 ns) are applied first, so the control region has the same selection as the signal region apart from the time.

Pitfall: the sign matters. `< -12.5` selects the early out-of-time side. `> -12.5` would select the in-time region, where the signal lives.

Note: the question text says "require cluster time < 12.5 ns", which conflicts with the instruction to use `-12.5`. The instruction and the code use `-12.5`, which is the only reading that gives an out-of-time region. This should be confirmed.

---

## Q: Correlation coefficient between |dphi| and N_hits

Question: "Compute the correlation coefficient for dphi and N_hits. The ABCD method works when the two chosen variables are independent for the background process."

```python
corr = scipy.stats.pearsonr(np.abs(deltaPhi), Nhits)[0]
print('Pearson correlation between |dphi| and Nhits: %.4f' % corr)
```

Explanation: `pearsonr` returns the coefficient and a p-value, and `[0]` keeps the coefficient. The absolute value is used because the ABCD regions are defined on |dphi|.

Expected result: a value close to zero, which supports the ABCD assumption. The value should be read from the notebook output.

Caveat: the Pearson coefficient only measures linear correlation. A small value is necessary but not sufficient, so the 2D plot and the scan below also test the assumption.

---

## Ex: 2D distribution of dphi and N_hits

Question: "Create a 2D histogram with axes named nhits and dphi, fill it with `Nhits` and `deltaPhi`, and plot it with `hist.plot2d(h, xaxis="nhits")`."

```python
h2 = hist.Hist("Events",
              hist.Bin("nhits", r"$N_{hits}$", 40, 0, 200),
              hist.Bin("dphi", r"$\Delta\phi(cls,MET)$", 64, 0, 3.2)
             )

h2.fill(nhits=Nhits, dphi=np.abs(deltaPhi))

hist.plot2d(h2, xaxis="nhits")
```

Explanation: the histogram has one axis per variable, filled with the per-cluster arrays. The axis is |dphi| because the ABCD regions use it. The binning puts the scan thresholds (N_hits 60 to 190, |dphi| 0.75) on bin edges.

What to look for: a tilt or diagonal band would mean the two variables are correlated, so the ABCD prediction would be biased. A flat pattern supports independence.

Difference from the question text: the text says to fill with `deltaPhi`. The code fills with `np.abs(deltaPhi)`, because the scan and the regions are defined on the absolute value.

---

## ABCD validation scan (N_hits or dphi)

Question: "Complete the number of events for bin B, C and D, and calculate the prediction and its uncertainty in bin A. Scan N_hits (scan_nhits = 1) or dphi (scan_nhits = 0)."

```python
        var = np.abs(deltaPhi)
        a = np.sum(np.logical_and(Nhits>=N_RECHIT_CUT, np.abs(deltaPhi)<DPHI_CUT))
        b = np.sum(np.logical_and(Nhits<N_RECHIT_CUT, var<DPHI_CUT))
        c = np.sum(np.logical_and(Nhits>=N_RECHIT_CUT, var>=DPHI_CUT))
        d = np.sum(np.logical_and(Nhits<N_RECHIT_CUT, var>=DPHI_CUT))

        pred = b*c/d
        unc_pred = pred*np.sqrt(1/b + 1/c + 1/d)
```

Explanation:
- A: high N_hits, small |dphi| (signal region; the number of observed events is already filled).
- B: low N_hits, small |dphi|.
- C: high N_hits, large |dphi|.
- D: low N_hits, large |dphi|.

If N_hits and |dphi| are independent for background, A/B = C/D, so the expected background in A is `pred = B*C/D`. The uncertainty treats each count as Poisson. For a product and a ratio the relative errors add in quadrature, so `unc_pred = pred * sqrt(1/B + 1/C + 1/D)`.

Pitfalls:
- Use the same |dphi| in all four regions. The `var` array does this.
- If D = 0, the prediction is infinite. This happens when the thresholds are too tight.
- With small counts the uncertainty is large, so agreement is not significant on its own.

---

## Discussion: ABCD method, closure and the meaning of the vetoes

ABCD method: the signal region A is defined by the two discriminating variables, N_hits (high) and |dphi| (small). Three control regions B, C and D are used to estimate the background in A, assuming the two variables are independent for background. Under that assumption the background in A is `B*C/D`, which is taken from data in B, C and D, so the estimate does not rely on simulation of the background.

Closure: the prediction is checked in regions where the signal is absent. Here the out-of-time region (`t < -12.5 ns`) is used. The scan over N_hits thresholds (60 to 190) and the |dphi| scan are performed, and the predicted value in A should agree with the observed value within the uncertainty for all thresholds. Agreement across the scan supports the independence assumption. A trend in the scan points to a correlation that would bias the signal region.

Meaning of the veto variables:
- Jet veto pT (GeV): the pT of the reconstructed jet matched to the cluster. A high value indicates a punch-through jet, so the cluster is background. Cut: pT < 10 GeV.
- Muon veto pT (GeV): the pT of the muon matched to the cluster. A high value indicates a muon bremsstrahlung shower. Cut: pT < 20 GeV.
- ME11/ME12 hits: hits in the first CSC stations. Background from jets and muon brems has these hits, signal is expected to have none.

Units: time in ns, pT in GeV, |dphi| in radians, N_hits dimensionless. The 2D histogram and the scan thresholds were chosen so that all thresholds sit on bin edges.

---

## Notes and open items

- The question file `_episodes/05-ABCD.md` is not in the repository. Questions are taken from the notebook.
- The OOT boundary (-12.5 ns) is taken from the notebook instruction. The episode 04 efficiency question uses an in-time window of -5 to 12.5 ns, so the two boundaries should be checked for consistency with the instructor.
- No code was run, so the printed values of the efficiencies, the correlation and the scan were not checked.
