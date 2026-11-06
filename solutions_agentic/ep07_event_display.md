# Episode 07: Event Display (optional)

Scope: `scripts/event_display.ipynb` (cells 0-7). The notebook was not executed.

## Inventory of blocks (episode 07-EventDisplay.md)

| # | Block | Line range | Status |
|---|---|---|---|
| 1 | Checklist: open `event_display.ipynb` | 19-24 | Covered by the notebook fill below |
| 2 | Shell step: `edmPickEvents.py` on the GEN-SIM-RECO dataset | 27-32 | Not a code blank; no fill needed |
| 3 | Copy-and-run `edmCopyPickMerge` command | 34-44 | Not a code blank; no fill needed |
| 4 | CernBox callout (Figure 5.1 image) | 46-52 | Not a question |
| 5 | cmsShow instructions and "Add Collection" note | 53-64 | Not a question |

Discussion and Question blocks: none in episode 07. This episode has no Discussion, Question, Exercise, Challenge or Solution block, so there are no questions to quote or answer.

Code blanks: none in `event_display.ipynb`. All cells are filled. One cell had a bug (below), which is now fixed.

## Code fill

### `event_display.ipynb`, cell 7 (write the event list)

Original problem: the cell computed `selection = np.logical_and(Nhits[k]>500, np.abs(deltaPhi[k])<0.75)` and printed its count, but then wrote `runNum[k][i]`, `lumiNum[k][i]`, `evtNum[k][i]` for the first ten events of the unmasked arrays. The events written therefore did not pass the cut the notebook prints.

New code:

```python
# select events
selection = np.logical_and(Nhits[k]>500, np.abs(deltaPhi[k])<0.75)
print(np.count_nonzero(selection))

# apply the selection mask before picking the first events
sel = np.asarray(selection, dtype=bool)
run_sel = np.asarray(runNum[k])[sel]
lumi_sel = np.asarray(lumiNum[k])[sel]
evt_sel = np.asarray(evtNum[k])[sel]

text_file = open(output_text_path + "/event_display.txt", "w")

# choose the first 10 selected events to view, don't change this number to more than 10,
# the next steps might take significantly longer
n_write = min(10, len(run_sel))
for i in range(n_write):
    print('{}:{}:{}'.format(run_sel[i], lumi_sel[i], evt_sel[i]))
    text_file.write('{}:{}:{} \n'.format(run_sel[i], lumi_sel[i], evt_sel[i]))
text_file.close()
```

Explanation:
- `edmPickEvents.py` expects one `run:lumi:event` per line. The file written by this cell is in that format.
- The mask is applied with `np.asarray(...)[sel]`, so the arrays are plain NumPy arrays regardless of whether uproot returned awkward arrays.
- `min(10, ...)` avoids an IndexError when fewer than ten events pass the cut.
- The cap of ten events is deliberate, because `edmPickEvents.py` and `edmCopyPickMerge` slow down quickly with more events (see the episode text).

Points for the instructor to decide (not changed):
- The threshold in the notebook is `Nhits > 500`. The episode says the events should "pass the signal region selections", which in episode 04 and in `create_datacard.ipynb` is `N_hits >= 130` and `|dphi| < 0.75`. If the display should show signal-region events, the `500` should be `130` (with `>=`). The 500 cut is kept because it is the value in the stripped notebook and selects large clusters, which are more visually interesting.

## Answers to discussion and question blocks

None. Episode 07 has no Discussion, Question, Exercise, Challenge or Solution blocks.

## Check status

- `python -I -m py_compile` passed on the new cell 7 source.
- The notebook was not executed. The `edmPickEvents.py` and cmsShow steps need a CMSSW environment and the CERNBox access described in the episode, and were not run.
