# Episode 04: Analysis strategy - filled code, explanations, and answers

Notebook: `analysis_strategy.ipynb` (exercises Ex0-Ex5). Inputs: episode `04-analysis_strategy.md`, stripped notebook.
Status: code filled from the lesson and notebook definitions. Numerical printouts were NOT produced here (notebooks not executed), so every numeric answer below says what to look for; the printed numbers from a run are the reference.

## Inventory (episode 04, line ranges in `04-analysis_strategy.md`)

| Block | Lines | Notebook cell(s) |
|---|---|---|
| Check your understanding (MET vs Higgs pT; trigger efficiency) | 58-66 | cells 6-7, 10-14 |
| Ex1 Discussion 4.1 (jet/muon veto) | 99-102 | cells 19-28 |
| Ex2 Discussion 4.2 (ME11/ME12 hits) | 118-122 | cells 29-34 |
| Ex3 Discussion 4.3 (cluster eta) | 138-142 | cells 35-39 |
| Ex4 Discussion 4.4 and 4.5 (cluster time) | 156-166 | cells 40-45 |
| Ex5 Question 4.1 + Solution 4.1 | 185-198 | cells 46-51 |
| Discussion 4.7 (N_hits and dPhi cuts) | 200-204 | cells 51-57 |

Note: there is no Discussion 4.6 in the episode text.

Notebook question blocks (all answered in the notebook as inserted markdown): cell 6 (Q: What do you observe?), cell 11 (signal efficiency), cells 22, 27, 33, 38, 43 (efficiencies), cell 45 (OOT region), cells 51, 56 (efficiencies).

## Units and conventions

- Cluster time in ns, eta dimensionless, dPhi in radians, ctau labels in mm (1000 mm = 1 m, 10000 mm = 10 m, 100000 mm = 100 m).
- Veto cuts are on pT in GeV (`< 10` jet, `< 20` muon).
- Empty events are filled with a sentinel (-1 or -999) so they fall into the underflow; `integral_low` therefore counts them in the numerator. This slightly inflates efficiencies for sparse samples. This is a known pitfall; keep it in mind when quoting numbers.
- Event-level masks are applied after the cluster mask: `X[sel_clusters][sel_ev]`.

## Filled code

### Ex0: MET and Higgs pT (cell 7)

```python
h= hist.Hist("Events",
              hist.Cat("var","Variable"),
              hist.Bin("pT", "MET/Higgs pT[GeV]", 50, 200, 1200)
             )

sample = "MC_40_1000"
# Fill the code
h.fill(var = "MET", pT = tree[sample].met )
h.fill(var = "Higgs pT", pT = tree[sample].gHiggsPt )

ax = hist.plot1d(h,density=True)
ax.legend(title = sample)
ax.set_ylabel("Fraction of events / bin width")
hep.cms.label(ax=ax,label="Preliminary")
```

### Ex0: Trigger efficiency (cells 10, 14)

```python
h_MET= hist.Hist("Events",
              hist.Cat("sample","sample"),
              hist.Bin("met", "MET[GeV]", 50, 0, 500)
             )

for m in mass:
    for ct in OLD_CTAU:
        # Fill the histogram h_MET with all the signals
        sample = 'MC_'+str(m)+'_'+str(ct)
        h_MET.fill(sample=sample, met=tree[sample].met)

ax = hist.plot1d(h_MET,density=True)
ax.set_ylabel("Fraction of events / bin width")
hep.cms.label(ax=ax,label="Preliminary")


for m in mass:
    for ct in OLD_CTAU:
        sample = 'MC_'+str(m)+'_'+str(ct)
        # Fill the code
        eff = integral(h_MET, sample, 200)
        print("sample = %s, eff(MET>200) = %.5f "%(sample, eff))
```

### Ex1: Jet veto pT (cells 20, 23)

```python
nBins    = 50
bin_lo   = 0
bin_high = 500

h = hist.Hist("Events",
              hist.Cat("sample","sample"),
              hist.Bin("pt", "Jet veto pT[GeV]", nBins, bin_lo, bin_high)
             )

for i, sample in enumerate(tree.keys()):

    #See if we are treating data
    isData = True if ("data" in sample) else False

    # build the selection arrays
    sel_clusters, sel_jet, sel_ev = select_events(tree[sample],isData)

    # get the Jet veto Pt
    VetoPt = tree[sample].cscRechitClusterJetVetoPt

    # Event and cluster selections (make sure we have a cluster [and matched to LLP in CSC] )
    # Cluster selection first (jagged), then the event selection on the outer axis
    VetoPt = VetoPt[sel_clusters][sel_ev]

    # Select the VetoPt of the largest cluster in the event
    # fill_none puts -1 for events without clusters that satisfy the cluster cuts
    VetoPt = ak.fill_none(ak.firsts(VetoPt,axis=1),-1 )

    h.fill(sample=sample,pt=VetoPt)


# Hint: use integral_low and integral

for i, sample in enumerate(tree.keys()):
    # Fill code here
    eff = integral_low(h, sample, 10)
    print("sample = %s, eff(VetoPT<10) = %.5f "%(sample, eff))
```

### Ex1: Muon veto pT (cells 25, 28)

```python
# Think about the bin ranges and number of bins you want to in the histogram before plotting
# 50 bins of 2 GeV: the 20 GeV threshold falls on a bin edge
nBins    = 50
bin_lo   = 0
bin_high = 100

h = hist.Hist("Events",
              hist.Cat("sample","sample"),
              hist.Bin("pt", "Muon veto pT[GeV]", nBins, bin_lo, bin_high)
             )


for i, sample in enumerate(tree.keys()):

    #See if we are treating data
    isData = True if ("data" in sample) else False

    # build the selection arrays
    sel_clusters, sel_jet, sel_ev = select_events(tree[sample],isData)

    # get the Muon veto Pt
    VetoPt = tree[sample].cscRechitClusterMuonVetoPt

    # Event and cluster selections (make sure we have a cluster [and matched to LLP in CSC] )
    VetoPt = VetoPt[sel_clusters][sel_ev]

    # Select the VetoPt of the largest cluster
    # fill_none puts -1 for events without clusters that satisfy the cluster cuts
    VetoPt = ak.fill_none(ak.firsts(VetoPt,axis=1),-1 )

    h.fill(sample=sample,pt=VetoPt)


# Hint: use integral_low and integral

for i, sample in enumerate(tree.keys()):
    # Fill code here
    eff = integral_low(h, sample, 20)
    print("sample = %s, eff(Muon VetoPT < 20) = %.5f "%(sample, eff))
```

### Ex2: ME11/ME12 hits (cells 31, 34)

```python
# Integer bins 0..10 so that the N=0 bin is [0,1)
nBins    = 10
bin_lo   = 0
bin_high = 10

h = hist.Hist("Events",
              hist.Cat("sample","sample"),
              hist.Bin("nME11_12", r"$N_{Me11/12 hits}$", nBins, bin_lo, bin_high)
             )


for i, sample in enumerate(tree.keys()):

    #See if we are treating data
    isData = True if ("data" in sample) else False

    # build the selection arrays
    sel_clusters, sel_jet, sel_ev = select_events(tree[sample],isData)

    # Compute nME11_12 (sum over the four ME11/ME12 chamber hit counts)
    nME11_12 = (tree[sample].cscRechitClusterNRechitChamberPlus11
              + tree[sample].cscRechitClusterNRechitChamberMinus11
              + tree[sample].cscRechitClusterNRechitChamberPlus12
              + tree[sample].cscRechitClusterNRechitChamberMinus12)

    # Event and cluster selections (make sure we have a cluster [and matched to LLP in CSC] )
    nME11_12 = nME11_12[sel_clusters][sel_ev]

    # Select the VetoPt of the largest cluster
    # fill_none puts -1 for events without clusters that satisfy the cluster cuts
    nME11_12 = ak.fill_none(ak.firsts(nME11_12,axis=1),-1 )

    h.fill(sample=sample,nME11_12=nME11_12)


# Hint: use integral_low and integral

for i, sample in enumerate(tree.keys()):
    # Fill code here
    eff = integral_low(h, sample, 1)
    print("sample = %s, eff(nME11_12 ==0) = %.5f "%(sample, eff))
```

### Ex3: Cluster eta (cells 36, 39)

```python
# 60 bins of width 0.1: the |eta| = 2 edges fall on bin edges
nBins    = 60
bin_lo   = -3
bin_high = 3

h = hist.Hist("Events",
              hist.Cat("sample","sample"),
              hist.Bin("eta", r"cluster $\eta$", nBins, bin_lo, bin_high)
             )


for i, sample in enumerate(tree.keys()):

    #See if we are treating data
    isData = True if ("data" in sample) else False

    # build the selection arrays
    sel_clusters, sel_jet, sel_ev = select_events(tree[sample],isData)

    # Get cluster Eta
    clsEta = tree[sample].cscRechitClusterEta

    # Event and cluster selections (make sure we have a cluster [and matched to LLP in CSC] )
    clsEta = clsEta[sel_clusters][sel_ev]

    # Select the eta of the largest cluster
    # fill_none puts -999 (not -1, which would fall inside |eta|<2) for events without clusters
    clsEta = ak.fill_none(ak.firsts(clsEta,axis=1),-999 )

    h.fill(sample=sample,eta=clsEta)


# Hint: use integral_low and integral
# |eta|<2 = (fraction below +2) - (fraction below -2); the underflow (-999) cancels in the difference

for i, sample in enumerate(tree.keys()):
    # Fill code here
    eff = integral_low(h, sample, 2) - integral_low(h, sample, -2)
    print(r'sample = %s, eff(abs($\eta$)<2) = %.5f '%(sample, eff))
```

### Ex4: Cluster time (cells 41, 44)

```python
nBins    = 100
bin_lo   = -25
bin_high = 25

h = hist.Hist("Events",
              hist.Cat("sample","sample"),
              hist.Bin("time", r"cluster time[ns]", nBins, bin_lo, bin_high)
             )


for i, sample in enumerate(tree.keys()):

    #See if we are treating data
    isData = True if ("data" in sample) else False

    # build the selection arrays
    sel_clusters, sel_jet, sel_ev = select_events(tree[sample],isData)

    # Get Variables
    clsTime  = tree[sample].cscRechitClusterTime
    jetVeto  = tree[sample].cscRechitClusterJetVetoPt
    muonVeto = tree[sample].cscRechitClusterMuonVetoPt
    clsEta   = tree[sample].cscRechitClusterEta

    # Apply the following selections:
    # - cluster and event selections
    # - JetVeto < 10 GeV
    # - MuonVeto < 20 GeV
    # - cluster eta < 2
    # Try running without Jet and muon Veto as well! (drop the jetVeto / muonVeto terms)
    sel_cls = sel_clusters & (jetVeto < 10) & (muonVeto < 20) & (np.abs(clsEta) < 2)
    clsTime = clsTime[sel_cls][sel_ev]

    # Select the time of the largest cluster
    # fill_none puts -999 for events without clusters that satisfy the cluster cuts
    clsTime = ak.fill_none(ak.firsts(clsTime,axis=1),-999 )

    h.fill(sample=sample,time=clsTime)


def integral_range(h,sample,thres_lo,thres_high):

    name = h.dense_axes()[0].name

    numer = h.integrate("sample",sample).integrate(name,slice(thres_lo,thres_high)).values()[()]
    denom = h.integrate("sample",sample).sum(name,overflow="over").values()[()]

    return numer/denom

for i, sample in enumerate(tree.keys()):
    # Fill code here
    eff = integral_range(h, sample, -5, 12.5)
    print(r'sample = %s, eff( -5 < t < 12.5 ns) = %.5f '%(sample, eff))
```

### Ex5: N_hits (cells 49, 52)

```python
nBins    = 30
bin_lo   = 0
bin_high = 300

h = hist.Hist("Events",
              hist.Cat("sample","sample"),
              hist.Bin("nhits", r"$N_{hits}$", nBins, bin_lo, bin_high)
             )


for i, sample in enumerate(tree.keys()):

    #See if we are treating data
    isData = True if ("data" in sample) else False

    # build the selection arrays
    sel_clusters, sel_jet, sel_ev = select_events(tree[sample],isData)

    # Get Variables
    nhits    = tree[sample].cscRechitClusterSize
    clsTime  = tree[sample].cscRechitClusterTime
    jetVeto  = tree[sample].cscRechitClusterJetVetoPt
    muonVeto = tree[sample].cscRechitClusterMuonVetoPt
    clsEta   = tree[sample].cscRechitClusterEta

    # Only plot Out-of-Time(OOT) region for data (t < -12.5 ns, early OOT)
    # For MC (signal) use the in-time signal window -5 < t < 12.5 ns
    if isData:
        sel_time = clsTime < -12.5
    else:
        sel_time = (clsTime > -5) & (clsTime < 12.5)

    # Apply the following selections:
    # - cluster and event selections
    # - JetVeto < 10 GeV
    # - MuonVeto < 20 GeV
    # - cluster eta < 2
    # - OOT selection for data; in-time selection for MC
    sel_cls = sel_clusters & (jetVeto < 10) & (muonVeto < 20) & (np.abs(clsEta) < 2) & sel_time
    nhits = nhits[sel_cls][sel_ev]

    # Select the value of the largest cluster
    # fill_none puts -1 for events without clusters that satisfy the cluster cuts
    nhits = ak.fill_none(ak.firsts(nhits,axis=1),-1 )

    h.fill(sample=sample,nhits=nhits)


for i, sample in enumerate(tree.keys()):
    # Fill code here
    eff = integral(h, sample, 130)
    print(r'sample = %s, eff($N_{hits}>130$) = %.5f '%(sample, eff))
```

### Ex5: |dPhi| (cells 54, 57)

```python
# 64 bins of width 0.05 on [0, 3.2]: the 0.75 threshold falls on a bin edge
nBins    = 64
bin_lo   = 0
bin_high = 3.2

h = hist.Hist("Events",
              hist.Cat("sample","sample"),
              hist.Bin("dphi", r"$\Delta\phi(cls,MET)$", nBins, bin_lo, bin_high)
             )


for i, sample in enumerate(tree.keys()):

    #See if we are treating data
    isData = True if ("data" in sample) else False

    # build the selection arrays
    sel_clusters, sel_jet, sel_ev = select_events(tree[sample],isData)

    # Get Variables
    dphi     = tree[sample].cscRechitClusterMet_dPhi
    clsTime  = tree[sample].cscRechitClusterTime
    jetVeto  = tree[sample].cscRechitClusterJetVetoPt
    muonVeto = tree[sample].cscRechitClusterMuonVetoPt
    clsEta   = tree[sample].cscRechitClusterEta

    # Only plot Out-of-Time(OOT) region for data, in-time window for MC
    if isData:
        sel_time = clsTime < -12.5
    else:
        sel_time = (clsTime > -5) & (clsTime < 12.5)

    # Apply the selections (cluster, event, vetoes, |eta|<2, OOT/in-time)
    sel_cls = sel_clusters & (jetVeto < 10) & (muonVeto < 20) & (np.abs(clsEta) < 2) & sel_time

    dphi = np.abs(dphi[sel_cls][sel_ev])

    # Select the value of the largest cluster
    # fill_none puts -999 for events without clusters that satisfy the cluster cuts
    dphi = ak.fill_none(ak.firsts(dphi,axis=1),-999 )

    h.fill(sample=sample,dphi=dphi)


for i, sample in enumerate(tree.keys()):
    # Fill code here
    eff = integral_low(h, sample, 0.75)
    print(r'sample = %s, eff($abs(\Delta\phi(cls,MET)<0.75)$) = %.5f '%(sample, eff))
```


## Explanations (what, why, intuition, pitfalls)

- **Trigger (Ex0).** No dedicated trigger exists for MDS in Run 2, so we use MET > 200 GeV. Signal passes only because the Higgs recoils against an ISR jet and the LLPs decay beyond the calorimeters (invisible in the calorimeter, so MET). The MET should track the generator Higgs pT at high values. Pitfall: MET is a sum over the event, so compare shapes, not event by event.
- **Jet and muon veto (Ex1).** Punch-through jets and muon brems showers from the main collision look like MDS clusters. Each cluster carries the pT of the jet/muon it was matched to; signal clusters have none, so vetoes at 10 GeV (jet) and 20 GeV (muon) remove most of the background with little signal loss. Pitfall: `integral_low` for these is an efficiency computed on the fill-none sentinel; read the efficiency for signal, not the total.
- **ME11/ME12 hits (Ex2).** These chambers sit behind little shielding, so showers that were not reconstructed as jets/muons still leave hits there. Requiring zero hits (`nME11_12 == 0`) is an active veto on those showers. Pitfall: the four chamber counts must be summed, and both plus/minus endcaps are needed.
- **Cluster eta (Ex3).** The muon reconstruction efficiency drops and pileup activity rises at large |eta|, so pileup clusters accumulate there. Cut at |eta| < 2. Pitfall: the sentinel for "no cluster" must not fall inside the window; I used -999 instead of the template's -1, because -1 lies inside |eta| < 2 and would bias the efficiency.
- **Cluster time (Ex4).** The LHC bunch spacing is 25 ns, so pileup clusters from neighbouring bunch crossings form peaks at roughly multiples of 25 ns, while signal and punch-through clusters peak in-time (t near 0). The signal window is -5 < t < 12.5 ns. The ECAL-free CSC timing resolution limits how sharp the peaks are.
- **Final variables (Ex5).** N_hits (cluster size) and |dPhi(cluster, MET)| are the two discriminants. Signal has large N_hits and dPhi peaked at 0; pileup clusters have small N_hits and flat dPhi. The dPhi and N_hits are independent for background, which is the ABCD premise in episode 05. Pitfall: plot the OOT region for data only, to avoid unblinding the signal region.

## Answers to the episode discussions and questions

Text is quoted verbatim from `04-analysis_strategy.md`.

**Check your understanding (lines 58-66):**
> For this part, open the notebook called `analysis_strategy.ipynb`
> You will
>  - confirm that pT of the Higgs matches the size of MET in the event.
>  - compute the trigger efficiency for signals
>
> The typical trigger efficiency for signals is around 1%

Answer: The MET and generator Higgs pT distributions have a similar high-pT shape, because the Higgs recoils against the ISR jet and the LLP decay products are invisible in the calorimeters. The trigger efficiency is the fraction of signal events with MET > 200 GeV; the lesson's ~1% is the typical expectation. The exact printed value from `analysis_strategy.ipynb` cells 14 (per sample) should be quoted; I did not run it.

**Discussion 4.1 (lines 99-101):**
> Does the shape make sense to you? Why does signal has small values and background has larger values? The analysis apply a cut of `cscRechitClusterJetVetoPt<10` and `cscRechitClusterMuonVetoPt<20`, do you agree with these selections?

Answer: Yes. Signal clusters from an LLP decaying in the muon system are not matched to a jet or a muon, so the veto pT is small (or -1 when no match). Background clusters from punch-through jets and brems showers are matched to energetic jets/muons, so their veto pT is large. The cuts at 10 GeV and 20 GeV keep signal essentially intact and remove the large-veto-pT tail. I agree with the selections; a looser cut would keep more background, a tighter one would start to cost signal.

**Discussion 4.2 (lines 118-120):**
> Does the shape make sense to you? Why does signal has small values and background has larger values? The analysis requires clusters to have no hits in ME11 and ME12, does you agree with the selection?

Answer: Yes. Showers from punch-through jets and brems deposit hits in ME11/ME12, which have little material in front of them. Signal clusters from the outer muon system rarely produce hits there. Requiring zero hits is an aggressive veto that removes the shower-like background that escaped the jet/muon vetoes. I agree, subject to checking the signal efficiency printed in cell 34.

**Discussion 4.3 (lines 138-140):**
> Does the shape make sense to you? The analysis requires clusters to have $ \| \eta \| < 2$, does you agree with the selection?

Answer: Yes. Pileup clusters are more numerous at large |eta| where muon coverage and reconstruction efficiency are lower, so the background rises there. The |eta| < 2 requirement removes this region with a modest signal loss, because the signal clusters are concentrated at central eta. I agree.

**Discussion 4.4 (lines 156-158):**
> Does the shape make sense to you? Why is data spaced at 25 ns?

Answer: The shape shows peaks at multiples of about 25 ns. The LHC bunch crossing period is 25 ns (40 MHz), so the cluster time is quantized by which bunch crossing the activity came from; the in-time peak is at t around 0 and out-of-time pileup peaks are at +-25 ns, -50 ns, etc. This is the signature that background clusters come from pileup interactions in other bunch crossings.

**Discussion 4.5 (lines 162-164):**
> Does the data distribution change with and without the vetoes applied? If it does, do you know why?

Answer: Yes, the in-time peak of the data distribution is reduced when the vetoes are applied. Punch-through jets and muon brems from the main collision are prompt (in-time), so they populate the t near 0 region; removing them reveals the out-of-time pileup structure. Run cell 41 with and without the jet/muon veto terms to see this.

**Question 4.1 (lines 185-187):**
> Why does the $\Delta\phi\text{(cluster, MET)}$ peak at 0 for signal, but flat for background distribution?

Solution (from `04-analysis_strategy.md` lines 190-197, reproduced in our words): For signal, the cluster lies along the LLP direction and the MET is along the Higgs direction, and the LLPs are back-to-back with the ISR jet, so the cluster and MET are aligned and dPhi peaks at 0. For background, clusters come from the underlying event while the MET comes from the primary hard interaction, so the two are uncorrelated and dPhi is flat. The flat dPhi also means dPhi is independent of N_hits for background, which is what the ABCD method requires.

**Discussion 4.7 (lines 200-202):**
> In the analysis, we apply a selection requiring the $N_{hits}>130$ and  $\| \Delta\phi\text{(cluster, MET)}\| < 0.75$. Do you agree with the selections?

Answer: Yes, with the caveat that the cut values should be checked against the efficiency printed in cells 52 and 57. N_hits > 130 selects the large, shower-like clusters that signal produces and removes small pileup clusters. |dPhi| < 0.75 retains the signal peak, which is broad because of the MET resolution and LLP direction, while removing the flat background. Both cuts are then the ABCD bin boundaries used in episode 05.
