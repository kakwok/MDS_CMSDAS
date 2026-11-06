# Episode 03 (MDS reconstruction): agentic solutions

Source: `scripts/MDS_reconstruction.ipynb`. The episode text `_episodes/03-MDS_reco.md` is not in this repository, so each question is quoted from the notebook prompt that asks for it.

---

## Exercise: compute cluster properties

> "Complete the function `computeCluster`, `computeStationProp` and `computeME11_12`" (with avgStation10 = sum_i(i x nhits_i) / sum_i(nhits_i), over stations with at least 10 hits)

```python
def computeCluster(rechits):
    cls = {}

    cls["Nhit"] = len(rechits)
    cls["Eta"] = sum(rechits.Eta)/cls["Nhit"]
    cls["Phi"] = sum(rechits.Phi)/cls["Nhit"]
    cls["X"] = sum(rechits.X)/cls["Nhit"]
    cls["Y"] = sum(rechits.Y)/cls["Nhit"]
    cls["Z"] = sum(rechits.Z)/cls["Nhit"]
    cls["Time"]= sum(rechits.Tpeak)/cls["Nhit"]
    ...

def computeStationProp(rechits):
    ...
        if nRechit_in_i_station >=10:
            nStation10   += 1
            avgStation10 += i_station*nRechit_in_i_station
            TotHitsStation10 += nRechit_in_i_station
    avgStation10 = avgStation10/TotHitsStation10
    return avgStation10,nStation10

def computeME11_12(rechits):
    # ME1/1 and ME1/2 are station +-1 with R < 465 cm (ME1/3 starts at R ~ 505 cm, see util.drawRZ)
    isME11_12 = np.logical_and(np.abs(rechits.Station)==1, (rechits.X**2+rechits.Y**2)**0.5 < 465)
    return sum(isME11_12)
```

Explanation:
- Each cluster property is the plain average over its rechits. Dividing by `Nhit` gives the (eta, phi) centroid used for the 0.4 cones in the event display, the global (x, y, z) centroid, and the time.
- `computeStationProp` loops over stations -4 to 4. Each station with at least 10 rechits adds 1 to `nStation10`, and adds `station * nhits` to the running sum. Dividing by the total hits in those stations gives the hit-weighted average station number.
- `computeME11_12` counts ME1/1 plus ME1/2 rechits. The ring is not in the ntuple branches used here, so the code selects station +-1 with R below 465 cm. The radius limits come from the ME1 boxes in `util.drawRZ`.
- Pitfalls: the `sum(...)/N` form avoids possible `np.mean` issues on awkward arrays. The plain phi average breaks for clusters that straddle plus or minus pi. The existing `avgStation10 = avgStation10/TotHitsStation10` line raises ZeroDivisionError when no station has at least 10 hits, because both values are Python integers. That line was not changed, so guard it if needed.
- Time uses `Tpeak` (strip-based peak time). `Twire` is the other option the episode mentions.

---

## Exercise: denominator selection

> "Compute the number of LLP that decays in the CSC first, then make a selection mask. How many events are there? This will be the denominator of your signal efficiency."

```python
denom_sel = np.logical_and(nLLPinCSC==1, signal.gLLP_csc)
```

Explanation:
- Event level: exactly one LLP decays in the CSC (`nLLPinCSC == 1`).
- Object level: the LLP itself is in the CSC (`signal.gLLP_csc`).
- Combining them broadcasts the event flag over the two LLPs, so `denom_sel` is a per-LLP mask that can index `gLLP_decay_vertex_z`.
- Pitfall: `gLLP_csc` alone would keep events where both LLPs decay in the CSC.

---

## Exercise: numerator selection

> "Select events with only 1 cluster AND in the denominator. Among those events, select the ones that there's one cluster matched to the LLP."

```python
numer_sel = np.logical_and(np.logical_and(nClsInCSC==1, nLLPinCSC==1), signal.cscRechitCluster_match_gLLP_csc)
```

Explanation:
- The numerator keeps the denominator events, requires exactly one matched CSC cluster (`nClsInCSC == 1`), and requires the matched cluster to be matched to a CSC-decaying LLP.
- The event-level conditions broadcast over the cluster axis, so `numer_sel` is a per-cluster mask.
- Pitfall: `nClsInCSC` counts matched clusters, not all clusters in the CSC. Dropping the `nLLPinCSC == 1` term changes the event sample.

---

## Exercise: inclusive efficiency

> "What is the overall efficiency of reconstructing an MDS cluster of 50 hits for this signal?"

```python
print("Events with 1 LLP decaying in CSC: ", ak.sum(denom_sel))
print("Events with 1 LLP decaying in CSC and a cluster matched to LLP: ", ak.sum(numer_sel))
print("Inclusive efficiency = %.2f"%(ak.sum(numer_sel)/ak.sum(denom_sel)))
```

Explanation:
- `ak.sum` of a boolean mask counts the True entries. Since each selected event has exactly one LLP, the denominator counts events.
- The efficiency is the numerator count divided by the denominator count for the 1 m (ctau = 1000 mm) signal sample with the 50-hit clustering.
- Pitfall: numerator and denominator must use the same event selection, otherwise the efficiency can exceed 1.

---

## Exercise: fill the histogram

> "Which variable should we fill-in in the function?"

```python
hz.fill(sample="denom", z = llp_z)
hz.fill(sample="numer", z = cls_z)
```

Explanation:
- `llp_z` is the decay z of the LLP in the denominator. `cls_z` is the same LLP decay z, but only for events where the LLP has a matched cluster.
- Events that fail a selection get zero, which lands in the underflow bin below 400 cm and is not drawn. This works only because the histogram range starts above zero.
- Pitfall: fill with the LLP decay position, not the cluster position, because the plot is efficiency as a function of where the LLP decayed.

---

## Exercise: 2D map of LLP and cluster positions

> "Make a 2D histogram to confirm your understanding!" (with the notebook asking for the LLP decay position in CSC and the matched cluster Z position, using the numerator selection)

```python
llp_z = ak.sum(abs(signal.gLLP_decay_vertex_z[denom_sel]),axis=1)
cls_z = ak.sum(abs(signal.cscRechitClusterZ[numer_sel]),axis=1)
hz.fill(z_llp=llp_z, z_cls = cls_z)
```

Explanation:
- Both arrays have one entry per event, with zeros where the event is not selected, so they can be filled event by event.
- Only events that pass both selections land in the drawn bins. Points on the diagonal mean the cluster sits at the same z as the LLP decay.
- Caveat: the cluster Z branch name `cscRechitClusterZ` is inferred from the `cscRechitClusterEta` and `cscRechitClusterPhi` naming used elsewhere, and was not verified against the ntuple.
- Pitfalls: take abs of the cluster z, since the histogram covers positive z only. The zero padding goes to underflow bins, which `hist.plot2d` does not draw.
