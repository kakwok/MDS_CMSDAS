# Episode 07: Event display (agentic solutions)

Scope: `scripts/event_display.ipynb`.

## Pick events from the signal region

> Save the run number, lumi section, and event number of the events that pass the signal region selection in a text file, then pick them with `edmPickEvents.py`.

### What was wrong

The skeleton computed `selection` but never used it. The loop wrote the first 10 events of the whole sample, not signal region events.

### Fix

```python
N_RECHIT_CUT = 130
DPHI_CUT = 0.75
selection = np.logical_and(Nhits[k]>=N_RECHIT_CUT, np.abs(deltaPhi[k])<DPHI_CUT)
print(np.count_nonzero(selection))

text_file = open(output_text_path + "/event_display.txt", "w")

for i in range(min(10, np.count_nonzero(selection))):
    print('{}:{}:{}'.format(runNum[k][selection][i], lumiNum[k][selection][i], evtNum[k][selection][i]))
    text_file.write('{}:{}:{} \n'.format(runNum[k][selection][i], lumiNum[k][selection][i], evtNum[k][selection][i]))
text_file.close()
```

Explanation for students:
- `selection` is the same signal region as bin A in the datacard (`Nhits >= 130` and `deltaPhi < 0.75`), so the events you view are the ones the analysis keeps.
- `runNum`, `lumiNum`, and `evtNum` are per-event arrays aligned with `Nhits` and `deltaPhi`, so the same boolean mask selects the matching entries.
- `min(10, ...)` guards against fewer than 10 selected events. The episode asks for at most 10 events because `edmPickEvents.py` gets slow with long lists.

Pitfalls:
- Indexing without `[selection]`, which silently shows the wrong events.
- Using `Nhits > 500` (the skeleton's earlier value) instead of the analysis cut, so the displayed events do not represent the selected signal region.
- Using a stale `k` from an earlier cell. It must be the signal key, `'MC_40_1000'` here.

### Effect on the episode commands

The `edmPickEvents.py` command and the `eventsToProcess=...` list in the episode were generated from the old, unfiltered list. After this fix, `event_display.txt` contains different events. Run `edmPickEvents.py` on the new `event_display.txt` and use the command it prints, rather than copying the list from the episode.

## Open the events with cmsShow

Not a coding step. After `edmCopyPickMerge` produces `pickevents.root`, copy it to a CERNBox directory with viewer access for `cms-vis-access`, open it at fireworks.cern.ch, and add the `csc2DRecHits` collection from the left panel. These collections are off by default.

## Physics sanity checks

- The event display is only a visual check. Use it to see whether the selected clusters look like the signal topology you expect (a displaced, large cluster in the CSC with a matching jet pattern).
- The selection uses reconstructed quantities from the ntuple, so the events shown are what the analysis would keep, not the generator-level LLP decay.
