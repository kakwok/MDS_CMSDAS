# Episode 05: ABCD background estimation - filled code, explanations, and answers

Notebook: `ABCD_validation.ipynb`. Inputs: episode `05-ABCD.md`, stripped notebook.
Status: code filled from the lesson and the notebook's own definitions. The notebook was NOT executed, so the agreement and correlation numbers must be read from a run.

## Inventory (episode 05)

| Block | Lines | Notebook cell(s) |
|---|---|---|
| ABCD method (equations and bin definitions) | 16-47 | cell 18 (text), cell 19 (code) |
| Open a notebook, OOT validation (checklist) | 56-59 | cells 4-5 |
| Discussion 5.1 (prediction vs observation) | 62-66 | cell 19 (scan) |

Notebook questions: cell 8 (Pearson correlation), cell 4 (OOT region text, `REPLACE_ME` at cell 5). Count of blanks filled: 4 code cells (5, 10, 14, 19).

## Units and conventions

- Time in ns (OOT cut at -12.5 ns; signal window -5 < t < 12.5 ns).
- N_hits are integer cluster sizes; dPhi in radians, |dPhi| used.
- Prediction: `A_pred = B * D / C`, with uncertainty propagated as relative Poisson errors `sigma_A/A = sqrt(1/B + 1/C + 1/D)`.

## Filled code

### Cell 5: OOT selection (REPLACE_ME -> -12.5)

```python
T = tree
        
########### SELECTION: CLUSTERS ############


sel_rechitcluster = np.abs(T['cscRechitClusterEta'].array()) < 2.0

me1112_veto = 0
sel_rechitcluster = np.logical_and(sel_rechitcluster, T['cscRechitClusterNRechitChamberPlus11'].array() <= me1112_veto)
sel_rechitcluster = np.logical_and(sel_rechitcluster, T['cscRechitClusterNRechitChamberPlus12'].array() <= me1112_veto)
sel_rechitcluster = np.logical_and(sel_rechitcluster, T['cscRechitClusterNRechitChamberMinus11'].array() <= me1112_veto)
sel_rechitcluster = np.logical_and(sel_rechitcluster, T['cscRechitClusterNRechitChamberMinus12'].array() <= me1112_veto)

sel_rechitcluster = np.logical_and(sel_rechitcluster, T['cscRechitCluster_match_MB1Seg_0p4'].array() <= me1112_veto)
sel_rechitcluster = np.logical_and(sel_rechitcluster, T['cscRechitCluster_match_RE12_0p4'].array() <= me1112_veto)
sel_rechitcluster = np.logical_and(sel_rechitcluster, T['cscRechitCluster_match_RB1_0p4'].array() <= me1112_veto)

sel_rechitcluster = np.logical_and( sel_rechitcluster, T['cscRechitClusterJetVetoPt'].array() < 10)
sel_rechitcluster = np.logical_and( sel_rechitcluster, T['cscRechitClusterMuonVetoPt'].array() < 20)
sel_rechitcluster = np.logical_and(sel_rechitcluster, T['cscRechitClusterTimeSpread'].array() <= 20)

# OOT selection
sel_rechitcluster = np.logical_and(sel_rechitcluster,  T['cscRechitClusterTime'].array() < -12.5)



########### SELECTION: JETS ############

sel_jet = np.logical_and(T['jetPt'].array() > 30, np.abs(T['jetEta'].array()) < 2.4 )
sel_jet = np.logical_and(sel_jet, T['jetTightPassId'].array())

########### SELECTION: EVENTS ############

sel_ev = T['METTrigger'].array()
sel_ev  = np.logical_and(sel_ev,np.sum(sel_rechitcluster,axis=1) >= 1)
sel_ev = np.logical_and(sel_ev ,T['met'].array() > 200)
sel_ev = np.logical_and(sel_ev , np.sum(sel_jet,axis=1)>=1)
sel_ev = np.logical_and(sel_ev,T['metFilters'].array())


##### cut based ID ####
cscRechitClusterNStation = T['cscRechitClusterNStation'].array()[sel_rechitcluster][sel_ev][:,0]
cscRechitClusterEta = T['cscRechitClusterEta'].array()[sel_rechitcluster][sel_ev][:,0]
cscRechitClusterPhi = T['cscRechitClusterPhi'].array()[sel_rechitcluster][sel_ev][:,0]
cscRechitClusterAvgStation = T['cscRechitClusterAvgStation'].array()[sel_rechitcluster][sel_ev][:,0]
cond2 = np.logical_and(np.abs(cscRechitClusterAvgStation)==2, np.abs(cscRechitClusterEta) < 1.6)
cond3 = np.logical_and(np.abs(cscRechitClusterAvgStation)==3, np.abs(cscRechitClusterEta) < 1.6)
cond4 = np.logical_and(np.abs(cscRechitClusterAvgStation)==4, np.abs(cscRechitClusterEta) < 1.8)
cond1 = np.logical_and(cscRechitClusterNStation==1, np.logical_or(np.logical_or(np.abs(cscRechitClusterAvgStation)==1, cond2), np.logical_or(cond3, cond4)))
cond2 = np.logical_and(cscRechitClusterNStation > 1, np.abs(cscRechitClusterEta) < 1.9)
bdt_sel = np.logical_or(np.logical_or(cond1, cond2), np.logical_or(cond3, cond4))
```

### Cell 10: Pearson correlation

```python
### print correlation between the two variables ####
# Pearson r and two-sided p-value for the background (out-of-time) clusters
rho, pval = scipy.stats.pearsonr(np.asarray(deltaPhi), np.asarray(Nhits))
print("Pearson rho(dPhi, Nhits) = %.4f  (p-value = %.3g, N = %d)" % (rho, pval, len(Nhits)))
```

### Cell 14: 2D Nhits vs dPhi

```python
# 2D histogram of Nhits vs |dPhi|, filled with the OOT background clusters
h2d = hist.Hist("Events",
                hist.Bin("nhits", r"$N_{hits}$", 50, 0, 200),
                hist.Bin("dphi", r"$\Delta\phi(cls,MET)$", 32, 0, 3.2))
h2d.fill(nhits=np.asarray(Nhits), dphi=np.abs(np.asarray(deltaPhi)))
hist.plot2d(h2d, xaxis="nhits")
```

### Cell 19: ABCD scan (b, c, d, pred, unc_pred filled)

```python
scan_nhits = 1

if scan_nhits:
    dphi_scan = [0.75]
    nhits_scan = np.arange(60,200,10)
else:
    dphi_scan = np.arange(9,22,1)/20
    nhits_scan = [130]


print('cut \t A \t B \t C \t D \t prediction \t uncertainty')
print()
for N_RECHIT_CUT in nhits_scan:
    for DPHI_CUT in dphi_scan:
        var = np.abs(deltaPhi)
        # A: high Nhits, small |dPhi|  (signal-like corner, observed)
        a = np.sum(np.logical_and(Nhits>=N_RECHIT_CUT, np.abs(deltaPhi)<DPHI_CUT))
        # B: high Nhits, large |dPhi|
        b = np.sum(np.logical_and(Nhits>=N_RECHIT_CUT, var>=DPHI_CUT))
        # C: low Nhits, large |dPhi|  (reference bin)
        c = np.sum(np.logical_and(Nhits<N_RECHIT_CUT, var>=DPHI_CUT))
        # D: low Nhits, small |dPhi|
        d = np.sum(np.logical_and(Nhits<N_RECHIT_CUT, var<DPHI_CUT))

        # ABCD prediction: A = B*D/C  (from N_A = c1*c2*C, N_B = c1*C, N_D = c2*C)
        pred = b*d/c if c>0 else np.nan
        # Poisson error propagation: (sigma_A/A)^2 = 1/B + 1/C + 1/D (relative, uncorrelated bins)
        unc_pred = pred*np.sqrt(1/b + 1/c + 1/d) if (b>0 and c>0 and d>0) else np.nan

        if scan_nhits:
            print(N_RECHIT_CUT, '\t',a,'\t',b,'\t',c,'\t',d,'\t', round(pred, 2), '\t',round( unc_pred, 2))
        else:
            print(DPHI_CUT, '\t',a,'\t',b,'\t',c,'\t',d,'\t', round(pred, 2), '\t', round( unc_pred, 2))
```


## Explanations (what, why, intuition, pitfalls)

- **OOT validation region.** The ABCD method needs Nhits and dPhi to be independent for background. The validation region is the early out-of-time (t < -12.5 ns) data: it is background-dominated (so it is signal depleted) and has the same pileup composition as the signal region. The check is that the ABCD prediction in this region matches what is observed.
- **Correlation (cell 10).** Pearson rho measures linear correlation. For background that is independent, rho is near 0. Pitfall: with large N, even a tiny rho can be "significant", so judge the magnitude as well as the p-value.
- **2D plot (cell 14).** A visual check for structure. A flat-looking 2D histogram supports independence. Pitfall: sparse high-N_hits bins look noisy; do not over-read them.
- **ABCD (cell 19).** Bins: A (high N_hits, small |dPhi|, signal-like), B (high N_hits, large |dPhi|), C (low N_hits, large |dPhi|, reference), D (low N_hits, small |dPhi|). From the lesson's model, N_A = c1 c2 Bkg_C, N_B = c1 Bkg_C, N_D = c2 Bkg_C, so A_pred = B*D/C. The c1 and c2 factors are the nuisance parameters for the N_hits and dPhi dependences. Pitfall: the notebook text gives `A = B x C / D`, which is inconsistent with that model; the code uses B*D/C, and this is flagged in the notebook.
- **Uncertainty.** Relative errors of B, C, D add in quadrature, and the prediction error is pred * sqrt(1/B + 1/C + 1/D). When a bin is empty the prediction is undefined; the code returns NaN in that case.

## Answers to the episode discussions and questions

Text quoted verbatim from `05-ABCD.md`.

**Discussion 5.1 (lines 62-64):**
> Does the prediction agree with the observation?

Answer (to be confirmed from the printed table of cell 19): The criterion is whether the observed A lies within the propagated uncertainty of B*D/C across the N_hits scan (dPhi < 0.75, N_hits from 60 to 190 in steps of 10). Agreement over the scan supports independence and background-only composition in the OOT region, which validates the method for use in the signal region. Expect larger uncertainties and possibly larger deviations at the tightest thresholds, where bins are small. A systematic deviation at high N_hits would indicate a residual correlation between the two variables for background, and the method would need a correction or a different choice of variables.

**Notebook question: Compute the correlation coefficient for dPhi and Nhits (cell 8):**
> The ABCD method works when the two chosen variables are independent of each others for *background* process.

Answer: Expect |rho| close to 0 for the OOT background sample; the printed value from cell 10 is the answer for your run. A rho that is clearly non-zero would be a warning that the independence assumption fails.
