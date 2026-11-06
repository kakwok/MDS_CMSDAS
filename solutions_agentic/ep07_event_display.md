# Episode 07 - Event display (optional): picking signal events and viewing them in cmsShow

This file documents the in-place solution of `scripts/event_display.ipynb` and the shell steps that follow it. The notebook itself is the "find the events" part; the three commands after it turn the list of event identifiers into a RECO-format ROOT file that the web-based `cmsShow` can open.

---

## 1. Ex: save run / lumi / event number of the events that pass the signal region

> This note book is very similar to the notebook in the previous episode, where we make the signal region selection. In addition, we save the run number, lumi section, and event number of the events that we want to view in a text file.

The path of the text file (already in the notebook):

```python
output_text_path = os.getcwd().replace('scripts', 'event_display')
if not os.path.isdir(output_text_path):os.makedirs(output_text_path)
print(output_text_path)
```

The final cell, as filled in:

```python
# select events
selection = np.logical_and(Nhits[k]>500, np.abs(deltaPhi[k])<0.75)
print(np.count_nonzero(selection))

text_file = open(output_text_path + "/event_display.txt", "w")

# choose the first 10 events to view, don't change this number to more than 10, 
# the next steps might take significantly longer
for i in np.nonzero(selection)[0][:10]:  # the first 10 events PASSING the signal region selection
    print('{}:{}:{}'.format(runNum[k][i], lumiNum[k][i], evtNum[k][i]))
    text_file.write('{}:{}:{} \n'.format(runNum[k][i], lumiNum[k][i], evtNum[k][i]))
text_file.close()
```

**Explanation.** The selection cell is the standard cluster + jet + event selection of the analysis (cluster `|eta| < 2.0`, no ME1/1 or ME1/2 hits, no matched muon-system segments, `jetVetoPt < 10`, `muonVetoPt < 20`, `timeSpread < 20`, `-5 < time < 12.5` ns; at least one `jetPt > 30` jet; `METTrigger`, `met > 200`, MET filters) plus the cut-based cluster ID, and it stores, per event, both the ABCD variables (`Nhits`, `deltaPhi`) and the identifiers `evtNum`, `runNum`, `lumiSec` of the DR event the ntuple row came from. On top of that the cell above applies `MC_condition == 2017`, which matters here and not in the analysis notebooks: the events we are going to pick live in a **2017** GEN-SIM-RECO dataset, so the run/lumi/event triplet must come from a 2017-mixing event, otherwise `edmPickEvents` will not find them.

The last cell then tightens the signal region to a very large cluster (`Nhits > 500`) pointing back at the MET (`|dphi| < 0.75`) - the "textbook" LLP decays that are worth looking at - and writes up to ten of them as `run:lumiSection:eventNumber` lines, the LUMI/EVENT list format that `edmPickEvents.py` reads. The one change with respect to the shipped cell is the loop header: `np.nonzero(selection)[0][:10]` returns the *indices of the events that pass the selection* (the first ten of them), whereas `range(10)` would have written the first ten events of the sample regardless of the tight requirement, i.e. the file and the printed `244` count would have described different events.

Pitfalls: the file is opened in `w` mode, so re-running the cell overwrites the previous list; keep the number of events at 10 or below (the dataset is scanned file by file over xrootd, so the pick time grows quickly); `Nhits`/`deltaPhi` and `runNum`/`lumiNum`/`evtNum` all have one entry per event (the leading selected cluster is taken with `[:,0]`), which is what makes the same index valid for both - if you ever change the notebook to loop over *all* clusters, that correspondence breaks.

---

## 2. Find the files that contain those events

```bash
edmPickEvents.py "/ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8/RunIIFall17DRPremix-PU2017_rp_94X_mc2017_realistic_v11-v1/GEN-SIM-RECO" event_display.txt
```

(the episode prints the dataset name with a typographic closing quote; use a plain `"`)

**Explanation.** `edmPickEvents.py` takes the dataset and the list of `run:lumi:event` triplets, queries DAS for the files of that dataset, and works out which file contains which requested lumi section. It does not copy anything itself - it prints an `edmCopyPickMerge` command that you then run. The dataset must be the same production the ntuples were made from (`RunIIFall17DRPremix ... 94X_mc2017_realistic_v11`), because the identifiers are only meaningful inside the run/lumi structure of that sample.

---

## 3. Copy the events into a local file

```bash
edmCopyPickMerge outputFile=pickevents.root \
  eventsToProcess=1:10122482,1:1055072,1:1055666,1:1127224,1:12138657,1:1441939,1:1441963,1:2514929,1:2896338,1:4255433 \
  inputFiles=/store/mc/RunIIFall17DRPremix/ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECO/PU2017_rp_94X_mc2017_realistic_v11-v1/130001/B2A708EC-0D4E-EA11-8B93-0025905C2CBC.root,/store/mc/RunIIFall17DRPremix/ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECO/PU2017_rp_94X_mc2017_realistic_v11-v1/130001/84696409-FC4D-EA11-A8CD-0CC47AFF0190.root,/store/mc/RunIIFall17DRPremix/ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECO/PU2017_rp_94X_mc2017_realistic_v11-v1/130000/349A8AB2-BF4D-EA11-9AFF-00259073E43C.root,/store/mc/RunIIFall17DRPremix/ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECO/PU2017_rp_94X_mc2017_realistic_v11-v1/130001/A4F52ADD-C94D-EA11-B339-7CD30AD091F0.root,/store/mc/RunIIFall17DRPremix/ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECO/PU2017_rp_94X_mc2017_realistic_v11-v1/130000/04DB0102-BB4D-EA11-A197-7CD30AD0A78C.root,/store/mc/RunIIFall17DRPremix/ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECO/PU2017_rp_94X_mc2017_realistic_v11-v1/130000/D6C63C14-484C-EA11-B04A-1866DA85DC8B.root,/store/mc/RunIIFall17DRPremix/ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECO/PU2017_rp_94X_mc2017_realistic_v11-v1/130001/4A6C137D-CA4D-EA11-B475-001E67E7195C.root,/store/mc/RunIIFall17DRPremix/ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RECO/PU2017_rp_94X_mc2017_realistic_v11-v1/130001/DAD03CA8-624E-EA11-A67F-0242AC1C0504.root
```

**Explanation.** This is the command printed in step 2, to be copy-pasted as is: `eventsToProcess` is `run:eventNumber` for each event to keep (note: **no** lumi section in this one, unlike the txt file) and `inputFiles` is the comma separated list of files that hold them. The result is `pickevents.root`, a small AOD/DR-format (`GEN-SIM-RECO`) file with exactly the events you selected. Your own run will print different files and possibly a different number of them, so always use the output of your own `edmPickEvents.py` rather than the text above.

---

## 4. Open the events with cmsShow

1. Create a new folder in your personal [CERNBox](https://cernbox.cern.ch), upload `pickevents.root`, and give **`cms-vis-access`** the *Viewer* permission on that folder (people -> `cms-vis-access`, otherwise the web viewer cannot read the file).
2. Open [fireworks.cern.ch](https://fireworks.cern.ch) (the web-based cmsShow), pick "open from CERNBox" and navigate to the file.
3. `csc2DRecHits` are **not** shown by default: use *Add Collection* on the left and add `csc2DRecHits` (and, if you want the full MDS picture, the other muon-system rechit collections). Then skim through the ten events.

**Explanation - why bother.** The whole analysis is built on a reconstructed "cluster" of CSC rechits, an object that exists only in the ntuple after the clustering code has run, so every selection studied so far is one level removed from the detector. The event display is the sanity check that cannot be done with ntuples: you want to see a compact cone/blob of 2D rechits in one or two stations of the muon system, with no charged track in the tracker pointing to it (that is what the `ME1/1`, `ME1/2`, `MB/RE/RB` match vetoes are rejecting), pointing back to the interaction point and to the MET direction, and sitting on top of the hadronic activity of the recoiling b-jets. If the surviving clusters instead look like single isolated hits, like a muon that crossed several stations in a straight line, or like noise spikes scattered over the wheel, then the "signal-like" yield is really detector background and the ABCD estimate of its rate would be the wrong tool. This is also the quickest way to understand why the `Nhits > 130` and time-window cuts do what they do.

Pitfalls: the web viewer shows the *reconstructed* collections of a `GEN-SIM-RECO` file, so tracker primitives and calorimeter towers are available but nothing from the MDS clustering itself (the clusters are ntuple-level objects, so what you see are the underlying `csc2DRecHits`); and the `CMS` label/overlay settings, plus switching to the 3D/Transverse view and hiding the muon "boxes", are usually needed before the shower is visible at all.
