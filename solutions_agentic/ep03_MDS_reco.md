# Episode 03 - MDS reconstruction (`scripts/MDS_reconstruction.ipynb`)

Solutions and explanations for the blanks in `MDS_reconstruction.ipynb`: Muon Detector Shower (MDS) rechits, the
DBScan clustering, the cluster properties (`computeCluster`, `computeStationProp`, `computeME11_12`), overlaying
signal on data, and the MDS reconstruction efficiency as a function of the LLP decay position.

---

## Ingredient of MDS: rechits

Rechits are loaded into one awkward-array object with `pack(events,"cscRechits")` (the `cscRechits` prefix is
stripped, so the fields are `Eta`, `Phi`, `X`, `Y`, `Z`, `Twire`, `Tpeak`, `Station`, `Chamber`, ...), giving
`d_rechits` for data and `s_rechits` for signal. A CSC rechit is the 2D hit made by crossing anode wires (one
coordinate + a time, `Twire`) with cathode strips (the orthogonal coordinate + a time, `Tpeak`) in the same gas
volume, which is why a single particle gives two coordinates.

### Discussion 3.1 - "Can you think of some of the properties of a rechit? Are there disadvantages/limitations for using rechit as the inputs of MDS?"

The properties that matter for MDS are **position** (`x,y,z` or equivalently `eta,phi`), **time** (`Twire`,
`Tpeak`) and the **location in the detector** (`Station` and `Chamber`, i.e. which endcap/station/substation the
chamber belongs to), plus simple quality/size information.

Limitation: the rechit reconstruction (local to the CSC, e.g. `CSCRecHit2D`) is **designed for a single muon** -
it associates wires and strips under the "one track" hypothesis. In a high-multiplicity shower that assumption
breaks down, so details of how the shower develops are lost (cluster shapes, per-muon separation, pulse-height
information). A dedicated ML reconstruction on the raw strips/wires could recover them. In practice the
limitation is minimal, because the most important MDS feature is the **multiplicity** (how many hits), and the
number of fired wires/strips is preserved even if the single-muon fit is imperfect.

---

## Clustering the rechits - DBScan and the cluster properties

DBScan is run on the `(eta, phi)` positions of each event's rechits with `minPts = 50`, `dR = 0.2` (50 is more
than twice the < 24 hits a single muon makes, so a 50-core cluster is not one muon), and `computeCluster` is
called on the constituent rechits of every cluster.

### Ex: `computeCluster` - position and time

```python
def computeCluster(rechits):
    cls = {}

    cls["Nhit"] = len(rechits)
    cls["Eta"] = np.mean(rechits.Eta)
    cls["Phi"] = np.mean(rechits.Phi)
    cls["X"] = np.mean(rechits.X)
    cls["Y"] = np.mean(rechits.Y)
    cls["Z"] = np.mean(rechits.Z)
    # cluster time: average of the wire time and the strip time of the rechits
    cls["Time"]= np.mean((rechits.Twire + rechits.Tpeak)/2.0)
    ...
```

**Explanation.** The input `rechits` is the *array of constituent rechits of one cluster in one event*, so every
property is just an average over that array (the definition given in the lesson: "cluster position = average of
input rechit positions", "cluster time = average of input rechit wire and strip times"). The mean position is a
good direction estimate for a compact shower blob; averaging the two independent times of every rechit first and
then averaging over rechits weights wire and strip measurements equally. A genuine shower is prompt, so its
average time lands on a bunch crossing, while accidentals give random times - which is what makes `Time` useful.

### Ex: `computeStationProp` - `nStation10` and `avgStation10`

```python
    for i_station in uniqueStations:          # [-4,-3,-2,-1,1,2,3,4]

        # total hits in the i-th station e.g. -4 
        nRechit_in_i_station = sum(rechits.Station==i_station)
        
        if nRechit_in_i_station >=10:
            # number of stations with >=10 hit
            nStation10   += 1
            # average of stations position with >=10 hit, weighted by the number of rechits in each station            
            avgStation10 += i_station*nRechit_in_i_station
            
            # total hits in station with at least 10 hits (i.e. sum of weights of avgStation10)
            TotHitsStation10 += nRechit_in_i_station
            
    avgStation10 = avgStation10/TotHitsStation10     # (already in the skeleton - see the bug note below)
```

**Explanation.** `nStation10` counts the stations that have at least 10 rechits of this cluster (`+= 1`), and
`avgStation10` accumulates the numerator of

$$\mathrm{avgStation10} = \frac{\sum_i i\times n_i}{\sum_i n_i},\qquad n_i = \text{rechits in station } i,\; n_i\ge 10,$$

so the division at the end completes the weighted average: `+= i_station*nRechit_in_i_station`. Requiring 10 hits
per station is what suppresses stations that were only touched by one or two unrelated muons. `Station` is
**signed** (positive/negative = +z/-z endcap), so `nStation10` runs 0-8 and `avgStation10` from -4 to +4: the
first tells you how much of the shower is spread over the muon system, the second tells you how deep it goes.

> **Bug (reported, left unchanged as instructed).** The skeleton line
> `avgStation10 = avgStation10/TotHitsStation10` (cell `12b8199a`, file line 440, inside `computeStationProp`)
> divides without a zero guard. If *no* station of a cluster has >= 10 rechits - possible even though the cluster
> has >= 50 rechits, e.g. hits spread over 6-8 stations with <= 9 each - then `TotHitsStation10 == 0` and Python
> raises `ZeroDivisionError`. The safe version would be
> `avgStation10 = avgStation10/TotHitsStation10 if TotHitsStation10 > 0 else 0`. It is outside the blanks to be
> filled here, so it was intentionally not modified (the reference solution has the same problem).

### Ex: `computeME11_12`

```python
def computeME11_12(rechits):
    # 'Chamber' encodes station and substation: +/-11 = ME1/1, +/-12 = ME1/2
    return ak.sum((abs(rechits.Chamber)==11) | (abs(rechits.Chamber)==12))
```

**Explanation.** ME1/1 and ME1/2 are *substations* of station 1, so `Station == 1` cannot select them - that is
what the `Chamber` branch is for: it is `station*10 + substation` with the sign of the endcap, so ME1/1 is
`+/-11` and ME1/2 is `+/-12`. Using `abs()` merges the two endcaps, and `ak.sum` over the boolean array counts
the rechits of the cluster that sit in those two innermost chambers. The analysis vetoes clusters with large
`nME11_12` because ME1/1-2 see a lot of punch-through from the calorimeters (a muon that wandered in from the
interaction region typically fires ME1/1 first).

### Ex: signal clusters and the overlay plots

```python
minPoint = 50
dR = 0.2
cls = DBScan(rechits,minPoint,dR)

# run the exact same clustering on the signal rechits (s_rechits, packed above)
s_cls = DBScan(s_rechits,minPoint,dR)
```

```python
samples=[
    {"cls":cls,"label":"data"},
    {"cls":s_cls,"label":"signal"},
    ...
]
```

**Explanation.** `s_rechits = pack(signal,"cscRechits")` is defined near the top of the notebook, so one extra
`DBScan` call with **identical** parameters produces `s_cls`, a cluster array with the same fields as `cls`;
appending it to `samples` is all that is needed because every plot cell loops over `samples` and fills the
`hist.Cat("sample")` category axis. Using the same `minPts`/`dR` for both is essential - otherwise the comparison
would contain a difference in reconstruction, not in physics. `ak.num(cls)` is an event-level quantity (number of
clusters per event) whereas `ak.flatten(cls.Nhit)` gives one entry per cluster. Switch `density=True` to compare
shapes rather than yields (the two samples have different numbers of events).

### Discussion 3.2 - "Which variables can be used to distinguish between signal and background?"

`N_rechits`, the cluster `Time` and `nME11_12` are the powerful ones:

- **Nhit**: an LLP shower produces hundreds of rechits, an accidental/cosmic single-muon cluster only tens, so
  the signal sits at much larger values. (This is the main variable, and the reason for `minPts = 50`.)
- **Time**: real showers are prompt, so `Time` peaks at a bunch crossing; background clusters assembled from
  unrelated hits give a wide, flat distribution. (The episode text defines the signal region as the early window
  `-12.5 < Time < -5` ns; note the notebooks' `create_datacard.ipynb` and one heading in
  `analysis_strategy.ipynb` instead use `-5 < t < 12.5` ns - see the mismatch reported in the ep04/ep06 write-ups.)
- **nME11_12**: background concentrates in the innermost chambers ME1/1-2, signal is spread over several
  stations, so requiring small `nME11_12` rejects a lot of background.
- `nStation10`/`avgStation10` and `Eta` are weak on their own but useful in combination (the "ClusterID").

Cuts on these variables give the signal/background separation used in the rest of the analysis.

---

## MDS reconstruction efficiency for signal

The efficiency is measured in the pre-clustered signal ntuple (`.../signal/nocuts/ggH_HToSSTobbbb_MH-125_MS-40_ctau-1000_...root`)
using only the fiducial sample "exactly one LLP decayed inside the CSC", so the result is a *pure reconstruction*
efficiency, not an acceptance times efficiency.

### Ex: Denominator selection

```python
nLLPinCSC = ak.sum(signal.gLLP_csc,axis=1)      # already in the skeleton

denom_sel = (nLLPinCSC==1) & (signal.gLLP_csc==1)
```

**Explanation.** `nLLPinCSC` sums the two LLP flags *within* each event, so `(nLLPinCSC==1)` is an **event-level**
selection: events with 0 CSC decays cannot be tested, events with 2 are ambiguous (we would not know which LLP a
cluster belongs to). The second factor, `signal.gLLP_csc==1`, is the **object-level** requirement "this LLP
decayed in the CSC", and it is what keeps the mask at LLP level, because it is applied as
`signal.gLLP_decay_vertex_z[denom_sel]` - we want to drop individual LLPs, not whole events. Awkward broadcasts
the event-level boolean over the two LLPs of the event, so the resulting `N * var` mask is `True` for exactly one
LLP in each selected event. `ak.sum(denom_sel)` is therefore the number of denominator events, and
`ak.sum(denom_sel)/len(signal)` is the CSC acceptance of the sample. `llp_z` is then the (absolute) decay z of
that LLP, one value per event.

### Ex: Numerator selection

```python
nClsInCSC = ak.sum(signal.cscRechitCluster_match_gLLP_csc,axis=1)   # already in the skeleton

numer_sel = (nClsInCSC==1) & (nLLPinCSC==1) & (signal.cscRechitCluster_match_gLLP_csc==1)
```

**Explanation.** Three conditions, as the cell comment spells out: exactly one cluster matched to a CSC-decaying
LLP (no ambiguity), the same event-level fiducial requirement as the denominator (`nLLPinCSC==1`, so the numerator
is a strict subset), and the cluster-level **match** flag itself, i.e. the cluster is within
$\Delta R(\text{cluster}, \text{LLP}) < 0.4$ of an LLP that decayed in the CSC. Note that the matching on the LLP
is already done for us in the ntuple (`cscRechitCluster_match_gLLP_*` branches), so we do not compute distances by
hand - we only require the match flag to be true and that there is exactly one such cluster. The mask is
cluster-level, which is why it can be applied to `signal.cscRechitCluster_match_gLLP_decay_z[numer_sel]` (the true
decay z of the matched LLP) and to `signal.cscRechitClusterZ[numer_sel]` (the measured cluster z).

### Ex: Inclusive efficiency

```python
print("Events with 1 LLP decaying in CSC: ", ak.sum(denom_sel))
print("Events with 1 LLP decaying in CSC and a cluster matched to LLP: ", ak.sum(numer_sel))
print("Inclusive efficiency = %.2f"%(ak.sum(numer_sel)/ak.sum(denom_sel)))
```

**Explanation.** Each mask holds at most one `True` per event, so summing all entries counts events, and the
ratio is the inclusive probability that an LLP decaying in the CSC is reconstructed as a DBScan cluster. It comes
out well below 1: a large fraction of the decays are missed - they do not make 50 hits close enough
together in eta-phi, e.g. because there is too little steel in front of the chamber where the LLP decayed.

### Ex: Fill the 1D histograms

```python
hz.fill(sample="denom", z = llp_z)
hz.fill(sample="numer", z = cls_z)
```

**Explanation.** `llp_z` and `cls_z` are event-level arrays (the `ak.sum(..., axis=1)` collapses the object
dimension), and both categories share the same axis (35 bins in 400-1075 cm = the length of the endcap muon
system, on `abs(z)`). Identical binning is what makes `hist.plotratio(hz.integrate("sample",'numer'),
hz.integrate("sample",'denom'), unc="clopper-pearson")` a legitimate binomial efficiency: each bin answers "of the
N decays at this z, how many were reconstructed?". Events that fail the selection carry z = 0, which is outside
the axis range and simply not filled.

### Ex: Plot the 2D histogram of LLP and cluster positions

```python
hz= hist.Hist("Events", 
              hist.Bin("z_llp", "LLP decay Z[cm]", 35, 400, 1075),
              hist.Bin("z_cls", "Cluster Z[cm]", 35, 400, 1075)              
             )
# Get the LLP decay position in CSC
llp_z = ak.sum(abs(signal.gLLP_decay_vertex_z[denom_sel]),axis=1)

# Plot the position of the matched cluster in z
# - Get the cluster Z position
# - Apply the numerator selection you computed above
cls_z = ak.sum(abs(signal.cscRechitClusterZ[numer_sel]),axis=1)

# Fill the 2D 
hz.fill(z_llp=llp_z, z_cls = cls_z)
```

**Explanation.** Now the axes are the **true** decay position and the **measured** cluster position, so this plot
checks where the reconstruction puts the shower. Only numerator events get a non-zero value in *both* coordinates
(denominator-only events have `z_cls = 0` and fall out of range), so the populated entries are exactly the
matched clusters. The band along the diagonal confirms that a cluster is reconstructed in the station where the
LLP decayed; off-diagonal entries are showers whose hits were picked up by a neighbouring (usually outer) chamber.

### Discussion 3.3 - "Why does the efficiency drop off at the two ends of the muon detectors? How does the efficiency vary between the muon stations?"

Two ingredients decide the efficiency: the LLP decay products need **some steel in front** to convert into a
shower of muons/hadrons, but the shower is only recorded in the **active gas** of the CSC chambers.

- **Front end (z < 568 cm, the `preME11` steel band of `drawCSCz`):** an LLP decaying there is still inside the
  calorimeters, in front of ME1/1. There is no CSC gas volume to record the decay and (for the part just in
  front of ME1/1) too little material converted the decay products yet, so few hits are made and no cluster
  passes 50 points.
- **Back end (z > 1075 cm, "Beyond CMS"):** the LLP has left the detector; there is no active chamber any more,
  so the efficiency falls to zero again.
- **Between stations:** the pattern follows the steel. The chambers sitting immediately **behind** a flux-return
  plate have the highest efficiency - ME1/2 behind ~39 cm of steel, ME2 behind ~65 cm, ME3 behind ~62 cm, ME4
  behind ~32 cm (the grey bands in `drawCSCz`/`drawCSCsteel`) - because the plate converts the decay products into
  a shower that then fires the gas gaps. Decays inside the gas volume *between* stations, and ME1/1 which has
  relatively little steel in front of it, produce smaller showers and are often lost. Hence the zig-zag
  (saw-tooth) structure of the efficiency vs z, and the 2D histogram above shows the reconstructed clusters
  concentrating in the stations just after the steel plates.

*Pitfalls to keep in mind:* the z of the *cluster* is not the z of the *decay* (the diagonal has a finite width),
and the efficiency is always relative to the fiducial selection - so this number cannot be compared to the
geometric acceptance of episode 02 without multiplying the two.
