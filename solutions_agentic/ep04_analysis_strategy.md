# Episode 04 — Analysis Strategy (`scripts/analysis_strategy.ipynb`)

This write-up goes through every blank that was filled in `scripts/analysis_strategy.ipynb`
(the trigger check, `Ex1`–`Ex5` and the efficiency questions) and answers the discussion
questions of episode *04-analysis_strategy* with the physics reasoning behind each selection.
Conventions used everywhere in the notebook: the ntuples are read with
`NanoEventsFactory.from_root(..., schemaclass=BaseSchema, treepath="MuonSystem")`, so every
branch is an awkward array (`tree[sample].branchName`), cluster branches are nested
(events → clusters) and event branches are flat; histograms are filled with `coffea.hist`
(`hist.Hist` / `hist.Cat` / `hist.Bin` / `hist.plot1d`) and the efficiencies are taken with the
provided helpers `integral`, `integral_low` and `integral_range`.

The full selection chain that the episode builds up is:

| step | requirement | kills |
|---|---|---|
| trigger | `MET > 200 GeV` (`met` + `METTrigger`) | everything not invisible enough |
| event | `metFilters`, ≥1 tight jet `pT>30`, `|eta|<2.4` | bad MET, non-ISR events |
| Ex1 | `cscRechitClusterJetVetoPt < 10`, `cscRechitClusterMuonVetoPt < 20` | punch-through jets, muon brem |
| Ex2 | `nME11_12 == 0` | punch-through/brem that were not reconstructed |
| Ex3 | `|cscRechitClusterEta| < 2` | forward pileup, poorly shielded region |
| Ex4 | `-12.5 < cluster time < -5` ns (SR), `time < -12.5` ns (OOT VR) | in-time background |
| Ex5 | `N_hits > 130`, `|Δφ(cluster,MET)| < 0.75` | ABCD discriminating variables |

---

## Trigger: plot the MET and the Higgs pT for signals

Question (notebook cell *Ex: Plot the MET and Higgs pT for signals*): "Focus on the region beyond
200 GeV — Q: What do you observe?"

```python
h= hist.Hist("Events",
             hist.Cat("var","Variable"),
             hist.Bin("pT", "MET/Higgs pT[GeV]", 50, 200, 1200)
            )

sample = "MC_40_1000"
# Fill the code
h.fill(var = "MET", pT = tree[sample].met )
h.fill(var = "Higgs pT", pT = tree[sample].gHiggsPt )
```

**Explanation.** `met` is the reconstructed MET and `gHiggsPt` the generator-level Higgs pT, both
flat per-event awkward arrays that `hist.fill` takes directly. The two shapes lie on top of each
other: in the higgs-portal benchmark the Higgs is produced back-to-back against one ISR jet, and if
both LLPs leave the calorimeters before decaying their momentum is invisible, so MET ≈ pT(H). This
is exactly the topology of Figure 4.2 of the episode and it is why a *generic* high-MET trigger can
be used to record an MDS signature. Pitfalls: the histogram is drawn with `density=True`, so only
shapes are comparable, and the axis starts at 200 GeV — everything softer goes into the underflow,
which is why the plot looks like a "trigger-side" sample.

---

## Trigger efficiency for signals

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
```

and

```python
for m in mass:
    for ct in OLD_CTAU:
        sample = 'MC_'+str(m)+'_'+str(ct)
        # Fill the code
        eff = integral(h_MET, sample, 200)
        print("sample = %s, eff(MET>200) = %.5f "%(sample, eff))
```

**Explanation.** `integral(h, sample, thres_lo)` integrates the dense axis from the threshold to
infinity (overflow included) and divides by the total, so with `200` it returns
`N(MET > 200)/N(all)` = the trigger efficiency of the signal. `200` GeV is a genuine bin edge of the
0–500 GeV, 50-bin axis (10 GeV wide), which is why the notebook prints
`h_MET.axis("met").edges()` first ("choose from the following bin edges"). The answer is of order
**1 %**, as the episode states. Physics: a 125 GeV Higgs decaying to two LLPs is normally soft, so
only the tail with a hard ISR jet recoil reaches 200 GeV of MET; longer-lived LLPs (larger ctau)
decay further outside the calorimeters, look more invisible, and therefore have a higher efficiency.
Pitfalls: (1) this is a *generator/reco-level* acceptance × trigger estimate on the `nocuts` samples
— no event selection is applied, on purpose; (2) because the samples are already
`..._weighted.root`, the ratio of weighted contents is the right number, do not re-weight by hand
here (that comes in the LLP lifetime episode).

---

## Ex0 — event/cluster selections (context)

`select_events()` (already filled in the notebook) returns `sel_clusters` (cluster level: at least
one cluster, and for MC one matched to a gen LLP that decayed in the CSC), `sel_jet` and `sel_ev`
(MET trigger, offline `met > 200`, `metFilters`, ≥1 tight jet, ≥1 selected cluster). Every exercise
below uses the same masking idiom that I filled in:

```python
VetoPt = VetoPt[sel_clusters][sel_ev]   # cluster-level mask, then event-level mask
```

**Explanation.** Applying `[sel_clusters]` first drops the clusters we do not want inside each event,
`[sel_ev]` then drops whole events; the comment in the cell about "the largest cluster" refers to
`ak.fill_none(ak.firsts(VetoPt,axis=1), -1)`, i.e. only the leading (largest) cluster of every event
is histogrammed, and events with no surviving cluster get a sentinel value (-1, or -999 for the time
plot) that falls in the underflow and therefore cancels out of the efficiency ratios.

---

## Ex1 — jet veto and muon veto pT

Jet veto (plot cells keep `50` bins from `0` to `500` GeV as given):

```python
    # get the Jet veto Pt
    VetoPt = tree[sample].cscRechitClusterJetVetoPt
    # Apply the cluster and event selections onto VetoPt
    VetoPt = VetoPt[sel_clusters][sel_ev]
```

Muon veto (the binning was a blank — `cscRechitClusterMuonVetoPt` is a pT in GeV, so I used the same
0–500 GeV range with 50 bins, which makes the 20 GeV working point a bin edge):

```python
nBins    = 50
bin_lo   = 0
bin_high = 500
...
    # get the Muon veto Pt
    VetoPt = tree[sample].cscRechitClusterMuonVetoPt
    # Apply the cluster and event selections onto VetoPt
    VetoPt = VetoPt[sel_clusters][sel_ev]
```

Efficiencies:

```python
    eff = integral_low(h, sample, 10)   # cscRechitClusterJetVetoPt < 10
    eff = integral_low(h, sample, 20)   # cscRechitClusterMuonVetoPt < 20
```

**Explanation.** `integral_low(h, sample, x)` = content below `x` / total content, i.e. the
efficiency of a `< x` cut; the thresholds 10 and 20 GeV are bin edges so the numbers are exact.

### Discussion 4.1

> "Does the shape make sense to you? Why does signal has small values and background has larger
> values? The analysis apply a cut of `cscRechitClusterJetVetoPt<10` and
> `cscRechitClusterMuonVetoPt<20`, do you agree with these selections?"

Yes, the shapes are exactly what the topology predicts. `cscRechitClusterJetVetoPt`
(`cscRechitClusterMuonVetoPt`) is the pT of the reconstructed jet (muon) that is matched to the
cluster. A signal cluster is a shower made by an LLP that decayed *outside* the calorimeters, so
there is nothing in the tracker or calorimeters pointing at it: no jet or muon is matched and the
variable sits at 0 (the small non-zero tail is accidental matching — a cluster that happens to lie
inside the cone of an unrelated ISR jet or of a real muon from the event). Punch-through jets and
muon bremsstrahlung, by contrast, *are* a jet and a muon: the cluster is the last visible part of a
genuine high-pT object, so the matched pT is large, up to several hundred GeV, which is the long tail
seen in data.

The cuts at 10/20 GeV are reasonable working points: they are placed where the background still has a
lot of events but the signal is already essentially at zero, so the loss of signal efficiency is tiny
while the bulk of the in-collision background is discarded. Too tight (e.g. 2–3 GeV) would start
throwing away signal from accidental matching and would also be dominated by jet energy/muon momentum
resolution effects; too loose (e.g. 50 GeV) would leave a large punch-through tail in the signal
region. The 20 GeV muon threshold is looser than the 10 GeV jet one because a muon that can radiate
enough bremsstrahlung to make a muon-system shower is typically harder, and because the muon pT scale
has larger tails than the jet energy scale, so a slightly higher threshold buys purity for very little
signal.

Pitfall worth remembering for the next episode: this veto only removes clusters matched to objects
from the *same* (hard-scatter) event. Pileup clusters are not matched to high-pT jets or muons, so the
vetoes barely touch them — which is why the out-of-time region stays background-enriched and why a
data-driven background estimate is needed at all.

---

## Ex2 — hits in ME11/ME12

```python
nBins    = 50
bin_lo   = 0
bin_high = 50

...
    # Compute nME11_12
    nME11_12 = (tree[sample].cscRechitClusterNRechitChamberPlus11  +
                tree[sample].cscRechitClusterNRechitChamberMinus11  +
                tree[sample].cscRechitClusterNRechitChamberPlus12  +
                tree[sample].cscRechitClusterNRechitChamberMinus12)
    # Apply the cluster and event selections
    nME11_12 = nME11_12[sel_clusters][sel_ev]
...
    eff = integral_low(h, sample, 1)   # nME11_12 == 0
```

**Explanation.** The four branches store the hits in ME1/+1, ME1/−1, ME1/+2, ME1/−2 separately, so
their sum is the total number of hits in the two least-shielded chambers. With 50 bins over 0–50 the
histogram has one bin per hit, hence "exactly zero hits" == "below 1 hit" and `integral_low(h, 1)` is
the efficiency of the requirement.

### Discussion 4.2

> "Does the shape make sense to you? Why does signal has small values and background has larger
> values? The analysis requires clusters to have no hits in ME11 and ME12, does you agree with the
> selection?"

Yes. ME1/1 and ME1/2 sit right behind the hadron calorimeter with little material in front, so
*anything that comes in from the inside* — a punch-through jet that escaped the calorimeters, a
bremsstrahlung photon/converting muon, and to some extent tracks from the primary interaction — crosses
one or both of them before it reaches ME2/ME3/ME4. Those background clusters therefore have a large
population at 1, 2, 3 … ME11/ME12 hits, even after the jet and muon vetoes (that is precisely why the
veto exists: badly reconstructed or sub-threshold objects leak through). A signal LLP, on the other
hand, decays *inside* the muon system, usually deep in it (the CSC stations at larger radius), so it
never deposits anything in the first station and the signal peaks at zero hits.

I agree with requiring `nME11_12 == 0`: it is a cheap, purely geometric ("active") veto that cleans
the leakage of the pT vetoes. The honest cost is signal efficiency: LLPs that decay early in the
endcap, or that decay in ME1/3 with hits spilling into ME1/2, are removed, and low-multiplicity
signal clusters are more likely to fail it, so the efficiency loss is lifetime dependent. Note also
that the strict `== 0` requirement is applied chamber-by-chamber in the analysis (and in
`ABCD_validation.ipynb`, which uses `<= 0` on each of the four branches plus the `MB1Seg/RE12/RB1`
match counters).

---

## Ex3 — cluster eta

```python
nBins    = 50
bin_lo   = 0
bin_high = 2.5

...
    # Get cluster Eta
    clsEta = abs(tree[sample].cscRechitClusterEta)
    # Apply the cluster and event selections
    clsEta = clsEta[sel_clusters][sel_ev]
...
    eff = integral_low(h, sample, 2)   # abs(eta) < 2
```

**Explanation.** I plot `|eta|` (the episode asks for `|eta|`, and the two muon endcaps are mirror
images), which makes the axis 0–2.5 with 0.05-wide bins so `|eta| < 2` is a bin edge and the efficiency
is a single `integral_low`.

### Discussion 4.3

> "Does the shape make sense to you? The analysis requires clusters to have |eta| < 2, do you agree
> with the selection?"

Yes. After the punch-through and brem vetoes the remaining background is pileup, and the pileup rate
per event rises with `|eta|`: there are more tracks/particles per unit eta in the forward region, the
muon system has less shielding and worse reconstruction/identification performance at high eta (so the
veto-based rejections become less efficient there). Data therefore rises steeply towards the edge of
the CSC acceptance (|eta| ≈ 2.4). Signal acceptance is comparatively flat — the LLP decay position is
set by the decay length and the geometry, not by the pileup density — and it falls only beyond the CSC
coverage. Requiring `|eta| < 2` keeps most of the signal while cutting the most pileup-rich region, and
it also keeps the cluster in the region where the ME11/ME12 shielding argument is strongest.

Pitfalls: cutting harder (e.g. `|eta| < 1.6`) starts to cost real signal acceptance for long-lived LLPs
that decay at large radius in the endcap disks, and it removes a large part of the OOT validation
region, which costs the background estimate statistics. The correct trade-off is where the signal
efficiency is still near its plateau — around 2.

---

## Ex4 — cluster time (after the vetoes)

```python
nBins    = 120
bin_lo   = -30
bin_high = 30

...
    # Get Variables
    clsTime  = tree[sample].cscRechitClusterTime
    jetVeto  = tree[sample].cscRechitClusterJetVetoPt
    muonVeto = tree[sample].cscRechitClusterMuonVetoPt
    clsEta   = abs(tree[sample].cscRechitClusterEta)

    # Apply the following selections:
    # - cluster and event selections
    # - JetVeto < 10 GeV
    # - MuonVeto < 20 GeV
    # - cluster eta < 2
    sel_cls = sel_clusters & sel_ev & (jetVeto < 10) & (muonVeto < 20) & (clsEta < 2)

    # Try running without Jet and muon Veto as well!
    clsTime = clsTime[sel_cls]
    # without the vetoes: clsTime = clsTime[sel_clusters & sel_ev]
...
    eff = integral_range(h, sample, -12.5, -5)   # signal region
```

**Explanation.** All four arrays are cluster level, so they can be AND-ed into one mask; `sel_ev` is
per event and awkward broadcasts the per-event boolean onto the cluster lists, which is why the
mixed-dimension `sel_clusters & sel_ev & (jetVeto < 10)` is legal. The binning is chosen so that
120 bins over −30…30 ns gives 0.5 ns bins and both −12.5 and −5 ns are bin edges (that is what makes
`integral_range` exact). `ak.fill_none(..., -999)` puts events without a surviving cluster in the
underflow.

Note: the question heading in the notebook reads "efficiencies for $-5<t<12.5$ ns", which is a typo in
the material (the sign of 12.5 is missing): the episode text says *"we apply a selection requiring the
cluster time to be between -5 ns and -12.5 ns as the signal region"*. I filled in the physics-correct
signal region `integral_range(h, sample, -12.5, -5)`; filling the printed range literally would
measure the in-time peak and would unblind the signal region for data.

### Discussion 4.4

> "Does the shape make sense to you? Why is data spaced at 25 ns?"

Yes. The LHC beams are bunched and consecutive bunch crossings are 25 ns apart, so pileup interactions
in the same fill happen at times 0, ±25, ±50 ns… relative to the hard scatter, and the CSC hit timing
is measured with respect to the hard-scatter (trigger) time. Clusters built from those other
interactions therefore pile up in peaks spaced by 25 ns instead of forming a smooth distribution — the
"comb" structure. (Prompt activity from the hard scatter itself is at t ≈ 0.) Signal is delayed rather
than quantised: an LLP travels a macroscopic distance before decaying, so its shower arrives late by an
exponentially distributed amount, and the signal distribution is broad and smooth.

The signal region −12.5 < t < −5 ns is the *gap* between the t = 0 bunch crossing and the previous one
at −25 ns: prompt particles from the hard scatter cannot be there (they arrive at ≈0 ns), and pileup
clusters from a neighbouring crossing can only leak into the edges of the gap (finite timing
resolution, wrong association of hits). A delayed LLP decay can be anywhere, including in this window,
so the window is nearly background free by construction and the background that remains is the smooth
pileup tail.

Pitfall: the very existence of the −25 ns peak is also a warning — if the window were defined sloppily
(e.g. −20 < t < −5) part of the out-of-time peak would enter the signal region.

### Discussion 4.5

> "Does the data distribution change with and without the vetoes applied? If it does, do you know why?"

Yes, and the way it changes is the strongest evidence that the analysis picture is right. Without the
jet/muon vetoes the distribution has a huge peak at t ≈ 0 ns (plus the 25 ns comb); with the vetoes ON
that in-time peak is largely gone, because it is made of punch-through jets and muon bremsstrahlung
from the main collision — objects that are always matched to a high-pT jet or muon. The out-of-time
peaks (|t| ≳ 10 ns) stay essentially unchanged: pileup clusters are not matched to any hard object, so
the vetoes do not see them.

That stability is also the *reason* the OOT region can be used as a validation region: if the vetoes
behaved differently in and out of time, the background composition of the OOT region would not be the
same as that of the signal region and the ABCD closure test would be meaningless. Run the cell with the
commented line (`clsTime[sel_clusters & sel_ev]`) to see the difference.

### Answer to "Can you suggest an out-of-time validation region?"

Take `cscRechitClusterTime < -12.5` ns, the *early out-of-time (OOT)* region: everything before the
signal region. It is signal depleted (an LLP decay would have to happen far earlier than the
hard scatter, which the exponential decay-time distribution practically never does — the signal time
distribution peaks at positive/delayed times), and it is populated by exactly the same out-of-time
pileup clusters, with the same N_hits and Δφ shapes, as the signal region, because the jet/muon/eta
vetoes act identically on both windows. A region further out (t < −30 ns) would be even cleaner but
statistically useless; a region sitting on the −25 ns peak would have a different background
composition.

---

## Ex5 — N_hits and Δφ(cluster, MET) in the OOT region

`N_hits`:

```python
nBins    = 50
bin_lo   = 0
bin_high = 500

...
    # Get Variables
    nhits    = tree[sample].cscRechitClusterSize
    clsTime  = tree[sample].cscRechitClusterTime
    jetVeto  = tree[sample].cscRechitClusterJetVetoPt
    muonVeto = tree[sample].cscRechitClusterMuonVetoPt
    clsEta   = abs(tree[sample].cscRechitClusterEta)

    # Only plot Out-of-Time(OOT) region for data
    if isData:
        sel_time = clsTime < -12.5
    else:
        sel_time = clsTime > -12.5

    sel_cls = sel_clusters & sel_ev & (jetVeto < 10) & (muonVeto < 20) & (clsEta < 2) & sel_time
    nhits = nhits[sel_cls]
...
    eff = integral(h, sample, 130)   # N_hits > 130
```

`|Δφ(cluster, MET)|` (same structure, `h = hist.Bin("dphi", ..., 60, 0, 3.0)`):

```python
    dphi     = abs(tree[sample].cscRechitClusterMet_dPhi)
...
    sel_cls = sel_clusters & sel_ev & (jetVeto < 10) & (muonVeto < 20) & (clsEta < 2) & sel_time
    dphi = dphi[sel_cls]
...
    eff = integral_low(h, sample, 0.75)   # |dphi| < 0.75
```

**Explanation.** Data is shown only in the OOT validation region (`t < -12.5` ns) so the signal region
is never unblinded; for the MC the complementary window (`t > -12.5` ns, i.e. signal region plus
in-time) is used so the signal shape is not artificially truncated by the same line of code. The
N_hits axis (0–500, 50 bins) and the |Δφ| axis (0–3.0, 60 bins) both have the working points 130 and
0.75 on bin edges; angles between 3.0 and π end up in the overflow, which `integral`/`integral_low`
count in the denominator (`sum(name, overflow="over")`), so the efficiencies stay correct. For data the
|Δφ| efficiency comes out close to the geometric value 0.75/π ≈ 0.24, a good sanity check that the
distribution really is flat.

Note that this exercise is labelled "Out-of-time Validation Region" because the ABCD background method
of the next episode is validated here, in the OOT region, not in the signal region.

### Question 4.1

> "Why does the Δφ(cluster, MET) peak at 0 for signal, but flat for background distribution?"

(Full text of this challenge question is in the episode page, lines 185–188 of
`/tmp/opencode/episodes/04-analysis_strategy.md`, with Solution 4.1 right below it at lines 190–198 —
the previous run reported that it could not be located; it is present in the gh-pages text.)

For signal, the cluster is the decay product of an LLP that comes out of the Higgs, so the cluster
direction follows the LLP boost direction, which is the Higgs direction. And the Higgs *is* the source
of the MET: it recoils against the ISR jet and nothing else in the event balances its momentum, so
`MET ≈ pT(H)` along the Higgs (the first plot of this notebook shows exactly this). Cluster direction ∥
Higgs ∥ MET ⇒ Δφ(cluster, MET) peaks at 0.

For background, the cluster is made of low-pT pileup particles from a *different* interaction, while
MET is computed from the primary event, so the two directions have nothing to do with each other and
Δφ is flat. A flat distribution in one variable means it carries no information about the other
variables that characterise the pileup cluster, in particular N_hits — and that independence of N_hits
and Δφ for background is precisely the assumption the ABCD method of episode 05 is built on.

Pitfall: "flat" is the expectation; a residual peak at Δφ ≈ 0 in the data would mean a remaining
non-pileup component (e.g. real MET-recoiling jets punching through), which is why the closure test is
done in the OOT region rather than assumed.

### Discussion 4.7

> "In the analysis, we apply a selection requiring the N_hits > 130 and
> |Δφ(cluster, MET)| < 0.75. Do you agree with the selections?"

Yes, with the caveat that these two cuts are not ordinary background-rejection cuts: they *define the
corners of the ABCD plane*, so they must be chosen for both purity and statistics.

* High multiplicity is the defining feature of the signature: an LLP decay produces a hadronic shower
  inside the chambers with hundreds of hits, whereas a pileup cluster is a random bunch of hits from
  low-pT particles and is almost always small. `N_hits > 130` therefore lands right on the signal peak
  and removes a huge amount of background (the data N_hits distribution in the OOT region falls
  steeply).
* `|Δφ| < 0.75` keeps the part of the plane where signal piles up (its Δφ peaks at 0) while remaining a
  compact box; for a flat background it simply keeps a fraction 0.75/π of the events.
* Both thresholds have to balance against statistics: the ABCD estimate needs a healthy number of
  events in *every* bin, in particular in the sidebands used for c1 and c2, and the OOT validation
  region has fewer events than the signal region by construction. Making the signal box tighter makes
  the closure test statistically meaningless.

The placement is optimised on things we are allowed to look at — the OOT data control sample and the
MC signal shape — and not on the signal region itself, which keeps the search unbiased. I agree with
the choice, with the note above that the values are a compromise between signal efficiency, background
rejection and ABCD statistical power.
