# Episode 03: MDS reconstruction (solutions, agentic)

**Discussion 3.1 (answers).** Question 1, verbatim: "Can you think of some of the properties of a rechit?" Question 2, verbatim: "Are there disadvantage/limitation for using rechit as the inputs of MDS?"

Q1: The properties most relevant for MDS are the position (x, y, z, eta, phi) and the time (Tpeak from the strips, Twire from the wires). The station is also useful, as are the chamber and ring that identify where the rechit sits.

Q2: Rechit reconstruction is designed for a single muon, so it can miss details of how a shower develops. A dedicated machine-learning reconstruction might recover those details. Multiplicity is the most important MDS feature, and rechits keep it, so the limitation is small.

## Challenge 1 - cluster properties (cell 24)

```python
def DBScan(rechits,nCore,eps):
    points=ak.zip([rechits.Eta,rechits.Phi],highlevel=False)

    all_clusters = []
    for x,rechit in zip(points,rechits):
        n_clusters_=0
        core_samples_mask=[]
        ## format into DBSCAN input format x=[ [x1,y1],[x2,y2] ... ]
        x=np.stack(ak.unzip(x),axis=1)
        i_clusters = []
        if len(x)>0:
            db = DBSCAN(eps=eps, min_samples=nCore).fit(x)
            core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
            core_samples_mask[db.core_sample_indices_] = True
            labels = db.labels_
            n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

            ## Loop over clusters in this event
            for iCs in range(n_clusters_):

                #Compute each cluster's properties:
                #    Input = rechits that belongs to this cluster (done via array masking)
                cls = computeCluster(rechit[labels==iCs])
                i_clusters.append(cls)
        all_clusters.append(i_clusters)
    return ak.Array(all_clusters)

# This function compute the cluster properties.
# The input `rechits` is an *array* of rechits used in the cluster
# This function is ran for each event
def computeCluster(rechits):
    cls = {}

    cls["Nhit"] = len(rechits)
    # position: average of the constituent rechit positions
    cls["Eta"] = np.mean(ak.to_numpy(rechits.Eta))
    cls["Phi"] = np.mean(ak.to_numpy(rechits.Phi))
    cls["X"] = np.mean(ak.to_numpy(rechits.X))
    cls["Y"] = np.mean(ak.to_numpy(rechits.Y))
    cls["Z"] = np.mean(ak.to_numpy(rechits.Z))
    # time: average of the wire time and the strip time of each constituent rechit
    cls["Time"] = np.mean(ak.to_numpy((rechits.Twire + rechits.Tpeak)/2.0))

    avgStation10,nStation10 = computeStationProp(rechits)
    nME11_12 = computeME11_12(rechits)

    cls["nStation10"] = nStation10
    cls['avgStation10'] = avgStation10
    cls["nME11_12"] = nME11_12

    return cls

# This function compute the cluster properties related to cluster stations
# The input `rechits` is an *array* of rechits used in the cluster in each event
# This function is ran for each event
def computeStationProp(rechits):
    # station can range from -4 to 4
    uniqueStations = np.array([-4,-3,-2,-1,1,2,3,4])
    nStation10   = 0
    avgStation10 = 0
    TotHitsStation10 = 0
    for i_station in uniqueStations:

        # total hits in the i-th station e.g. -4
        nRechit_in_i_station = sum(rechits.Station==i_station)

        #
        if nRechit_in_i_station >=10:
            # number of stations with >=10 hit
            nStation10   += 1
            # average of stations position with >=10 hit, weighted by the number of rechits in each station
            avgStation10 += i_station*nRechit_in_i_station

            # total hits in station with at least 10 hits (i.e. sum of weights of avgStation10)
            TotHitsStation10 += nRechit_in_i_station

    # guard: a cluster with no station having >=10 hits would divide by zero
    avgStation10 = avgStation10/TotHitsStation10 if TotHitsStation10 > 0 else -999.0
    return avgStation10,nStation10

# This function compute the total number of rechits in ME11 and ME12
# The input `rechits` is an *array* of rechits used in the cluster
# This function is ran for each event
def computeME11_12(rechits):
    # ME11 = station 1 ring 1, ME12 = station 1 ring 2
    # NOTE: the 'Ring' field name is assumed; confirm with rechits.fields (cell 15)
    is_ME1 = np.abs(rechits.Station) == 1
    is_ring12 = (rechits.Ring == 1) | (rechits.Ring == 2)
    return int(ak.sum(is_ME1 & is_ring12))
```

Each cluster property is a direct transcription of the definitions in the lesson:
- Nhit = number of constituent rechits.
- Eta, Phi, X, Y, Z = plain averages of the rechit positions.
- Time = average of (Twire + Tpeak)/2 over the rechits.
- nStation10 = number of stations with at least 10 rechits.
- avgStation10 = sum(i * nhits_i)/sum(nhits_i) over the stations with at least 10 rechits.
- nME11_12 = number of rechits in station 1 with ring 1 or 2.

Pitfalls: (1) Phi is a naive average and does not handle the wrap at +/-pi, which is fine for the clusters here but worth knowing. (2) avgStation10 divides by zero if no station has 10 rechits; a guard now returns -999 for that case, which the histograms drop as out of range. (3) The ring field name is an assumption; confirm it in rechits.fields.

## Challenge 2 - signal clusters (new cell after 25)

```python
s_cls = DBScan(s_rechits,minPoint,dR)
```

**Challenge 2 (explanation).** The signal rechits (s_rechits, built in cell 13) are clustered with the same DBSCAN parameters (minPts = 50, dR = 0.2) and stored in s_cls, so the signal and data clusters go through identical code.

## Challenge 2 - add signal to samples (cell 31)

```python
samples=[
    {"cls":cls,"label":"data"},

    # you may add more lines like this one to overlay them on plots
    # e.g. {"cls":s_cls,"label":"signal"},
    # where `s_cls` is the clusters from signals
    {"cls":s_cls,"label":"signal"},
]
density=False
```

**Challenge 2 (explanation, plots).** Adding the signal line to `samples` overlays signal (s_cls) and data (cls) in every plot below. Compare N_rechits, Time, and nME11_12 in particular: these are the variables that separate signal from background in the next discussion. Reading the plotting code, each cell fills a Hist with a `sample` category axis and draws the normalised or raw distributions side by side.

**Discussion 3.2 (answer).** Question, verbatim: "Which variables can be used to distinguish between signal and background?"

The distributions of N_rechits, Time, and ME11_12 differ strongly between signal and background. Signal clusters are made of many rechits, have a Time that is spread in a distinct way, and carry more ME1/1 and ME1/2 rechits. Cuts on these variables separate signal from background. nStation10 and Eta are weaker on their own, which is why the bonus cell combines nStation10 with eta.

## Ex: Denominator selection (cell 52)

```python
denom_sel = (nLLPinCSC == 1) & signal.gLLP_csc
```

nLLPinCSC is an event-level count, so (nLLPinCSC == 1) is one boolean per event. Multiplying by signal.gLLP_csc, which is one boolean per LLP, keeps only the LLPs that decayed in the CSC, in events where exactly one LLP did. Pitfall: an event-level array must not be used as the object mask on its own; the & broadcasts it onto the LLP level.

**Denominator answer.** The count of events with one LLP in the CSC is printed in cell 53. It is the denominator of the efficiency. Divide it by the total number of generated events to get the CSC acceptance.

## Ex: Numerator selection (cell 57)

```python
numer_sel = (nClsInCSC == 1) & (nLLPinCSC == 1) & signal.cscRechitCluster_match_gLLP_csc
```

Three conditions: exactly one cluster matched to a CSC LLP (event-level), exactly one LLP in the CSC (event-level), and the cluster itself is matched to a CSC LLP (object-level). The first two are broadcast onto the cluster level by the &. nClsInCSC counts clusters matched to a CSC LLP, as defined in cell 56.

## Ex: Inclusive efficiency (cell 60)

```python
print("Events with 1 LLP decaying in CSC: ", ak.sum(denom_sel))
print("Events with 1 LLP decaying in CSC and a cluster matched to LLP: ", ak.sum(numer_sel))
print("Inclusive efficiency = %.2f"%(ak.sum(numer_sel)/ak.sum(denom_sel)))
```

**Inclusive efficiency (answer).** The notebook question is "What is the overall efficiency of reconstructing an MDS cluster of 50 hits for this signal?"

The efficiency is the number of events with one CSC LLP and a matched cluster, divided by the number of events with one CSC LLP. Both counts are printed by cell 60. Because each selected event has exactly one LLP and at most one matched cluster, ak.sum counts events directly. I did not execute the notebook, so the printed value is not quoted here.

## Fill the histogram (cell 62)

```python
#Create a 2D histogram with a "sample" axis and an axis
hz= hist.Hist("Events",
              hist.Cat("sample","sample"),
              hist.Bin("z", "LLP decay Z[cm]", 35, 400, 1075)
             )

#Fill the LLP decay position in CSC in the denominator bin
hz.fill(sample="denom", z = llp_z[llp_z>0])
#Fill the LLP decay position in CSC matched to a cluster in the numerator bin
hz.fill(sample="numer", z = cls_z[cls_z>0])
```

**Notebook question (answer).** Question, verbatim: "Which varaible should we fill-in in the function?"

For the denominator, fill the LLP decay Z of the events with one LLP in the CSC (llp_z from cell 54). For the numerator, fill the LLP decay Z of the matched clusters (cls_z from cell 58). Both are in cm. The filter `>0` drops the zero entries that appear for events outside the selection, which would otherwise land in the underflow bin.

**Discussion 3.3 (answer, part 1 and 2).** Questions, verbatim: "Why does the efficiency drops off at the two ends of the muon detectors?" and "How does the efficiency varies between the muon stations? Do you understand why?"

The efficiency falls at the two ends for different reasons. At large z, beyond the last CSC station (ME4), there is no detector behind the LLP decay, so a shower is not seen. At small z, before ME1/1, the decay happens in the thick steel of the yoke, where the shower can be absorbed before it produces enough rechits to form a cluster.

Between the stations, the efficiency is modulated by the steel. Decays just in front of a steel layer develop a shower that reaches the chambers, while decays inside a gas gap or right after a steel layer produce less. Read the station-by-station pattern off the plot and compare it with the steel positions in drawCSCz. I did not run the notebook, so the exact values are not quoted.

## Ex: 2D histogram of LLP and cluster z (cell 69)

```python
hz= hist.Hist("Events",
              hist.Bin("z_llp", "LLP decay Z[cm]", 35, 400, 1075),
              hist.Bin("z_cls", "Cluster Z[cm]", 35, 400, 1075)
             )
# Event mask: events that have a matched cluster (numerator events)
ev_numer = ak.sum(numer_sel, axis=1) > 0

# Get the LLP decay position in CSC, for the numerator events
llp_z = ak.sum(abs(signal.gLLP_decay_vertex_z[denom_sel]), axis=1)[ev_numer]

# Plot the position of the matched cluster in z
# - Get the cluster Z position
# - Apply the numerator selection you computed above
# NOTE: the cluster Z branch name is an assumption; confirm it with signal.fields (cell 48)
cls_z = ak.sum(abs(signal.cscRechitClusterZ[numer_sel]), axis=1)[ev_numer]

# Fill the 2D
hz.fill(z_llp=llp_z, z_cls = cls_z)
```

The 2D histogram puts the LLP decay Z on one axis and the matched cluster Z on the other. Entries lie near the diagonal when the cluster sits at the decay point. Off-diagonal spread shows that the shower develops over several cm of steel, so the cluster Z is shifted downstream of the decay. ev_numer keeps the two arrays the same length.

**Discussion 3.3 (answer, part 3).** Question, verbatim: "Make a 2D histogram to confirm your understanding!"

The 2D histogram shows the cluster Z against the LLP decay Z for matched events. A band along the diagonal means clusters form near the decay. Its width and the shift of the band show how far the shower develops before the cluster is formed. The pattern supports the efficiency interpretation above: decays at positions that do not give a developed shower in the chambers appear as a gap in the band.

