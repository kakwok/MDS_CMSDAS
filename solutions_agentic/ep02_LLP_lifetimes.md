# Episode 02 (LLP lifetimes): agentic solutions

Source: `scripts/LLP_lifetimes.ipynb`. The episode text `_episodes/02-LLP_intro.md` is not in this repository, so each question is quoted from the notebook prompt that asks for it.

Units: the tree branches are in cm, the file names and `ctau_list` are in mm, and the plots are in cm or m. Each solution states which unit it uses.

---

## Ex0: LLP proper decay length

> "Fill in the righthand side of each equation below with the 4 branches mentioned to calculate the proper decay length"

```python
gLLP_ctau = {}
for k, T in tree.items():
    gLLP_decay_vertex = np.sqrt(T['gLLP_decay_vertex_x'].array()**2 + T['gLLP_decay_vertex_y'].array()**2 + T['gLLP_decay_vertex_z'].array()**2)
    gLLP_gamma = 1.0/np.sqrt(1.0 - T['gLLP_beta'].array()**2)
    gLLP_ctau[k] = gLLP_decay_vertex/(T['gLLP_beta'].array()*gLLP_gamma)
```

Explanation:
- The lab-frame decay length is the distance from the primary vertex (the origin) to the decay vertex, `sqrt(x^2 + y^2 + z^2)`, in cm.
- `d = beta * gamma * ctau`, so `ctau = d / (beta * gamma)`. `gamma = 1/sqrt(1 - beta^2)` comes from `gLLP_beta`.
- The result is the proper decay length (`c * tau`) in cm.
- Pitfalls: build gamma from the velocity (`gLLP_beta`), not from a momentum. The LLP branches have shape (events, 2), so keep the operations element-wise. The decay length assumes the LLP is produced at the origin, which holds for these generator-level samples. Remember that 1000 mm in a file name is 100 cm in the branches.

---

## Ex1: plot and fit the proper decay length

> "replace the `REPLACE_ME` variable below to fill in the correct variable to extract the fitted ctau by using the variable `fitted_slope`" (the notebook also says the file names are in mm and the branches are in cm)

```python
    fit = h[k].GetFunction('expo')
    fitted_slope = fit.GetParameter(1) # expo is exp(p0 + p1*x), so p1 = -1/ctau, in cm
    print("sample: {}, fitted ctau [mm]: {}".format(k, -10.0/fitted_slope))
```

Explanation:
- The histogram of `gLLP_ctau` is exponential, so it is a straight line on the log scale.
- ROOT's `expo` model is `exp(p0 + p1 * x)`, so `p1 = -1/ctau`. The answer is `-1/p1` in cm.
- The label asks for mm, so the code multiplies by 10 (`-10.0/p1`). A 1000 mm sample should print about 1000.
- Pitfalls: the slope is negative, so keep the minus sign. The factor is 10, not 1/10. The bins are 100 cm wide, which is coarse for the 1 m sample, so expect rough agreement there.

---

## Ex2: reweighting function

> "You will implement this equation in the function `weight_calc` below"

```python
def weight_calc(llp_ct, new_ctau, old_ctau):
    denominator = np.exp(-1.0*llp_ct/old_ctau)/old_ctau**2
    numerator = np.exp(-1.0*llp_ct/new_ctau)/new_ctau**2
    weight = numerator/denominator
    
    return weight
```

Explanation:
- The old per-event probability is `denominator`, and the new one is `numerator`. Their ratio is `w = [exp(-t/tau_new)/tau_new^2] / [exp(-t/tau_old)/tau_old^2]`. This is the `(tau_old/tau_new)^2 exp[t (1/tau_old - 1/tau_new)]` form in the notebook.
- `llp_ct` is `t1 + t2`. The two LLPs decay independently, so the joint probability is a product of exponentials and depends only on the sum.
- Pitfalls: all three inputs must use the same units (cm). Reweighting to a much longer lifetime than the source makes a few events dominate. `weight_calc` in `lib/helper.py` implements the same formula and can be used as a cross-check.

---

## Reweight the 10 m sample to 1 m

> "You will need to fill in the arguments for the `weight_calc()` function with `llp_ct` and the two lifetimes. Be sure to convert everything to cm."

```python
        llp_ct = np.sum(gLLP_ctau[k], axis = 1)
        
        # Fill in the arguments in weight_calc()
        weight = weight_calc(llp_ct, 100.0, 1000.0) # new = 1 m = 100 cm, old = 10 m = 1000 cm
                
        for j in range(Nevents_todraw): h[k].Fill(gLLP_ctau[k][j][0], weight[j])
```

Explanation:
- The target is 1 m (`tau_new = 100` cm). The source is the 10 m sample (`tau_old = 1000` cm).
- The argument order is `(llp_ct, new_ctau, old_ctau)`. Swapping new and old gives the inverse weight.
- Each event carries one weight, so the histogram is filled with the first LLP's ctau and `weight[j]` for event `j`. The fitted lifetime should be close to 1000 mm, matching the native 1 m sample.

---

## Ex3: geometric decay region in the muon system

> "Replace the `REPLACE_ME` in the cell below." (the cell asks for the endcap region based on Fig. 4.1.1 of the Muon TDR; the probability should be about 15-20% for the 1 m and 10 m samples)

```python
    sel_ms = np.abs(T['gLLP_eta'].array()) < 2.4
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_z'].array()) < 1075)
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_z'].array()) > 400)
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_r'].array()) < 700)
    
    print(k, np.count_nonzero(np.sum(sel_ms, axis = 1))/len(sel_ms))
```

Explanation:
- A LLP counts as decaying in the endcap muon system if its decay vertex is inside `|eta| < 2.4` (muon coverage), `400 < |z| < 1075` cm (the endcap station range; 1075 cm is where the drawings end), and `r < 700` cm (outer radius of the muon system).
- `np.abs` on z includes both endcaps. `np.sum(sel_ms, axis=1)` counts how many of the two LLPs pass.
- Physics: the decay probability depends on the lab-frame decay length, about `beta * gamma * c * tau`, so the acceptance rises with the lifetime.
- Caveat: these boundaries are a choice from the geometry and were not checked against the figure. A rough estimate with a 1 m sample (about 100 cm) gives only a few percent of LLPs decaying beyond 400 cm, so the 15-20% target may refer to the 10 m sample. Compare the printed numbers with that expectation.

---

## Acceptance versus proper decay length

> "You will fill in the definition of `sel_ms`, which is what you defined in the cell above and the arguments for `weight_calc()` to reweight to intermediate lifetimes in the `ctau_list`"

```python
    sel_ms = np.abs(T['gLLP_eta'].array()) < 2.4
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_z'].array()) < 1075)
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_z'].array()) > 400)
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_r'].array()) < 700)
    llp_ct = np.sum(gLLP_ctau['MC_40_{}'.format(ct_source)], axis = 1)
    weight = weight_calc(llp_ct, ct/10.0, ct_source/10.0) # mm -> cm
```

Explanation:
- `sel_ms` repeats the Ex3 cuts on the sample named by `ct_source`.
- `ct_source` is the shortest generated lifetime that is at least as long as the target. Reweighting to a shorter lifetime keeps the weights bounded. The last point (500 m) has no longer sample, so it uses the 100 m sample and has larger weight fluctuations.
- The target and source lifetimes are converted from mm to cm before calling `weight_calc`.
- Acceptance is `sum(weight over selected events) / sum(weight over all events)`. Normalising by the number of events instead would be wrong.
