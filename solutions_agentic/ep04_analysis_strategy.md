# Episode 04: Analysis strategy (agentic solution)

Notebook: `scripts/analysis_strategy.ipynb`. The question file `_episodes/04-analysis_strategy.md` is not in this repository or on the `2025` branch, so the questions below are quoted from the notebook markdown cells. The code was not executed, so no numbers are quoted. The answers give the expected behaviour, which should be checked by running the notebook.

Common conventions used in every fill:
- Per-cluster variables are masked with the cluster mask first and the event mask second, e.g. `var[sel_clusters][sel_ev]`.
- `ak.firsts(..., axis=1)` keeps the first cluster per event, and `ak.fill_none(..., value)` marks events without a surviving cluster.
- Units: MET, pT and veto pT in GeV; cluster time in ns; eta and |dphi| dimensionless and radians.

---

## Ex: MET and Higgs pT for signal (Q: What do you observe?)

Question: "Focus on the region beyond 200 GeV. Q: What do you observe?"

```python
h.fill(var = "MET", pT = tree[sample].met )
h.fill(var = "Higgs pT", pT = tree[sample].gHiggsPt )
```

Explanation: both variables are per event, so they fill directly. The `var` category keeps the two curves separate, and `density=True` normalizes them.

Expected observation: above 200 GeV both distributions fall steeply and have similar shapes. The signal MET is the recoil of the Higgs against the boosted ISR jet, so a large MET requires a large Higgs pT. This is why the MET trigger is a reasonable proxy for the signal topology.

Pitfall: events below 200 GeV are in the underflow bin of this histogram, so they are not visible.

---

## Q: Signal efficiency for MET > 200 GeV

Question: "What is the signal efficiencies for MET>200 GeV?"

```python
eff = integral(h_MET, sample, 200)
```

Explanation: `integral` returns the number of events with MET at or above the threshold, including overflow, divided by all events of the sample. The 200 GeV threshold is a bin edge of the 10 GeV-wide histogram, so the result is exact.

Expected observation: the efficiency is a modest fraction of the signal, because the large-MET tail is populated only by high-pT Higgs events. Read the printed value from the notebook.

---

## Ex1: Jet veto and muon veto (Q: efficiency of jet veto pT < 10 GeV)

Question: "What are the efficiencies Jet veto pT < 10 GeV?" (muon veto analogous, below)

```python
VetoPt = VetoPt[sel_clusters][sel_ev]
...
eff = integral_low(h, sample, 10)
```

Explanation: the jet veto pT is the pT of the reconstructed jet matched to the cluster. Keeping clusters with a small value removes punch-through jet background, because a punch-through jet deposits energy that the calorimeter does not fully absorb, so the matched jet has high pT.

Expected observation: signal clusters are not matched to jets, so most of them have a veto pT near zero. The efficiency for pT < 10 GeV should therefore be high for signal and lower for the background-enriched data sample.

Pitfall: the first cluster is chosen by `ak.firsts`, so the order of clusters in the ntuple matters.

## Q: Muon veto plot and efficiency (muon veto pT < 20 GeV)

Question: "What are the efficiencies Muon veto pT < 20 GeV?"

```python
nBins    = 50
bin_lo   = 0
bin_high = 500
...
VetoPt = VetoPt[sel_clusters][sel_ev]
...
eff = integral_low(h, sample, 20)
```

Explanation: the muon veto pT is the pT of the muon matched to the cluster. A high-pT muon can produce a bremsstrahlung shower in the muon system, so clusters matched to a high-pT muon are removed by a veto below 20 GeV.

Binning: 10 GeV bins from 0 to 500 GeV put the 20 GeV threshold on an edge.

Expected observation: the muon veto distribution is strongly peaked near zero, as the log scale of the plot suggests. The 20 GeV requirement removes only a small fraction of the signal.

## Ex2: ME11/ME12 hits (Q: efficiency of nME11_12 == 0)

Question: "What are the efficiencies for nME11_12 ==0?"

```python
nME11_12 = (tree[sample].cscRechitClusterNRechitChamberPlus11 + tree[sample].cscRechitClusterNRechitChamberMinus11
            + tree[sample].cscRechitClusterNRechitChamberPlus12 + tree[sample].cscRechitClusterNRechitChamberMinus12)
nME11_12 = nME11_12[sel_clusters][sel_ev]
...
eff = integral_low(h, sample, 1)
```

Explanation: the four chamber-hit counts are summed per cluster. Bins are one hit wide with integer edges, so bin 0 holds exactly the clusters with no ME11/ME12 hits. `integral_low(h, sample, 1)` integrates up to the edge at 1.

Physics: background from jets and muon brems has hits in the first CSC stations, which have little shielding in front of them. Signal clusters from LLP decays in the CSC are expected to have none.

Pitfall: integer counts need integer bin edges. Non-integer edges would split the zero count between bins.

## Ex3: Cluster eta (Q: efficiency for abs(eta) < 2)

Question: "What are the efficiencies for abs(eta)<2 ?"

```python
nBins    = 50
bin_lo   = -2.5
bin_high = 2.5
...
clsEta = tree[sample].cscRechitClusterEta
clsEta = clsEta[sel_clusters][sel_ev]
clsEta = ak.fill_none(ak.firsts(clsEta,axis=1),-999 )
...
eff = integral_low(h, sample, 2) - integral_low(h, sample, -2)
```

Explanation: `integral_low` gives the fraction below a threshold. The window `-2 < eta < 2` is the difference of two such fractions, with the same denominator. The underflow entries cancel in the difference.

Change from the template: the template fills events without a cluster with -1. That value lies inside the eta range and would show up as a spike at eta = -1 and bias the efficiency. The fill value is -999, which is in underflow.

Expected observation: the high-|eta| region carries more background (lower muon reconstruction efficiency and more pileup), so the |eta| < 2 cut keeps most signal and removes a larger fraction of background.

## Ex4: Cluster time after the vetoes (Q: efficiency for -5 < t < 12.5 ns)

Question: "What are the efficiencies for -5<t<12.5 ns ?"

```python
clsTime  = tree[sample].cscRechitClusterTime
jetVeto  = tree[sample].cscRechitClusterJetVetoPt
muonVeto = tree[sample].cscRechitClusterMuonVetoPt
clsEta   = tree[sample].cscRechitClusterEta
...
sel_cls = sel_clusters & (jetVeto < 10) & (muonVeto < 20) & (abs(clsEta) < 2)
clsTime = clsTime[sel_cls][sel_ev]
clsTime = ak.fill_none(ak.firsts(clsTime,axis=1),-999 )
...
eff = integral_range(h, sample, -5, 12.5)
```

Explanation: all vetoes are combined into one cluster mask so that every variable refers to the same clusters. The time bins are 0.5 ns wide from -20 to 40 ns, so the edges -5 and 12.5 ns are exact.

Expected observation: signal clusters from the in-time bunch crossing lie inside the window, while the pileup background is spread over a wider range. The efficiency of the in-time window is therefore high for signal and lower for background.

## Q: Out-of-time validation region (discussion)

Question: "Can you suggest an out-of-time validation region based on the above plot? A good validation region needs to have: no (or small) signal contamination (signal depleted), and a similar background composition as the signal-rich region."

Answer: choose clusters with time below the in-time window, `t < -5 ns` at the lower edge of the in-time window seen in the time plot, or further out, as done in `ABCD_validation.ipynb` with `t < -12.5 ns`.

Reasoning:
- Signal depleted: the signal is produced in time, so the signal clusters peak inside the window. Going further out in time removes the tail of signal clusters near the boundary.
- Similar background: the out-of-time clusters come from the same detector and pass the same vetoes, so their composition (pileup-dominated, low pT) matches the in-time background.

Caveat for the notebooks: the in-time window in the efficiency question is -5 to 12.5 ns, while the ABCD notebook uses -12.5 ns as the OOT boundary. This solution uses -12.5 ns in both places so the two notebooks agree. The boundary should be confirmed by the instructor.

## Ex5: Number of hits N_hits (Q: efficiency for N_hits > 130)

Question: "What are the efficiencies for N_hits > 130 ?"

```python
nhits    = tree[sample].cscRechitClusterSize
...
sel_cls = sel_clusters & sel_time & (jetVeto < 10) & (muonVeto < 20) & (abs(clsEta) < 2)
nhits = nhits[sel_cls][sel_ev]
...
eff = integral(h, sample, 130)
```

Explanation: `integral(h, sample, 130)` counts clusters with N_hits at or above 130. The question writes `>130`. For integer counts the ABCD notebook uses `>=`, so the same convention is kept. Bins are 5 hits wide, so 130 is an edge.

Blinding: data are filled only out of time (`t < -12.5 ns`). MC is filled in time (`t > -12.5 ns`).

Expected observation: signal clusters have a larger number of hits than the background, so the efficiency for N_hits >= 130 is higher for signal than for the data control sample.

## Ex5: |dphi(cluster, MET)| (Q: efficiency for |dphi| < 0.75)

Question: "What are the efficiencies for abs(dphi(cls,MET)<0.75) ?"

```python
dphi     = np.abs(tree[sample].cscRechitClusterMet_dPhi)
...
dphi = dphi[sel_cls][sel_ev]
...
eff = integral_low(h, sample, 0.75)
```

Explanation: the absolute azimuthal separation between the cluster and the MET is used. The histogram has 0.05 rad bins from 0 to 3.2, so 0.75 is a bin edge.

Expected observation: signal clusters are expected to lie close in phi to the MET direction, so the efficiency of |dphi| < 0.75 is high for signal and lower for background.

Pitfall: the variable is in radians.

---

## Question 4.1 (challenge)

Not resolved. The challenge question is not in the notebook and the `_episodes` file is not in the repository, so its text could not be read. Nothing is written for it.
