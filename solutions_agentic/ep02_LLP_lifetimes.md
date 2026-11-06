# Episode 02 - LLP lifetimes (`scripts/LLP_lifetimes.ipynb`)

Solutions and explanations for every blank in the `LLP_lifetimes.ipynb` notebook: calculating the LLP proper
decay length, fitting it, reweighting the lifetime, and computing the geometric acceptance of the endcap muon
system as a function of the mean proper decay length.

General convention for this episode (and the biggest source of errors): **the ROOT file names are in mm**
(`ctau-1000` = 1000 mm = 1 m) while **all tree branches are in cm** (`gLLP_decay_vertex_z`, `gLLP_ctau`, ...).

---

## Ex0 - "Fill in the righthand side of each equation below with the 4 branches mentioned to calculate the proper decay length"

```python
gLLP_ctau = {}
for k, T in tree.items():
    # lab-frame decay length: distance from IP (0,0,0) to decay vertex, in cm
    gLLP_decay_vertex = np.sqrt(T['gLLP_decay_vertex_x'].array()**2 + T['gLLP_decay_vertex_y'].array()**2 + T['gLLP_decay_vertex_z'].array()**2)
    gLLP_beta = T['gLLP_beta'].array()
    # gamma = 1/sqrt(1-beta^2)
    gLLP_gamma = 1.0/np.sqrt(1.0 - gLLP_beta**2)
    # proper decay length: ctau = d/(beta*gamma), in cm
    gLLP_ctau[k] = gLLP_decay_vertex/(gLLP_gamma*gLLP_beta)
```

**Explanation.** Particle decay is a Poisson process, so in the particle's own rest frame the survival
probability is exponential, $N(t) = e^{-t/\tau}$. The simulation, however, only records where the LLP decayed in
the **lab**, and the lab decay length is stretched by time dilation:

$$d = \beta\gamma\, c\tau \quad\Longrightarrow\quad c\tau = \frac{d}{\beta\gamma}.$$

1. `gLLP_decay_vertex` $=\sqrt{x^2+y^2+z^2}$ is the distance from the interaction point (which is at the origin
   of CMS coordinates) to the decay vertex, in cm.
2. `gLLP_beta` is $\beta=v/c$ and `gLLP_gamma = 1/sqrt(1-beta**2)` is the Lorentz factor.
3. Dividing by $\beta\gamma$ "de-boosts" the measurement to the LLP rest frame. The result, $c\tau$, is a
   **Lorentz invariant** property of the particle: unlike $d$, it does not depend on how the LLP happened to be
   boosted in a given event, and therefore it should be exponentially distributed with the mean that was set at
   generation time.

*Pitfalls:* forgetting $\beta$ (i.e. dividing only by $\gamma$), or mixing units - the output of this cell is in
**cm**, not mm. Each event has two LLPs, so the branches have shape `(events, 2)`.

---

## Ex1 - "Replace the `REPLACE_ME` variable below to fill in the correct variable to extract the fitted ctau"

```python
    h[k].Fit('expo')
    fit = h[k].GetFunction('expo')
    print("sample: {}, fitted ctau [mm]: {}".format(k, -10./fit.GetParameter(1)))
```

**Explanation.** ROOT's `'expo'` fit function is $\ln N = p_0 + p_1 x$, so the slope is $p_1 = -1/\tau$ in the
units of the fitted variable. Here the fitted variable is the histogram axis, i.e. `gLLP_ctau` in **cm**:

$$c\tau[\mathrm{cm}] = \frac{-1}{p_1},\qquad c\tau[\mathrm{mm}] = \frac{-10}{p_1}.$$

`fit.GetParameter(1)` returns $p_1$, hence `-10./fit.GetParameter(1)`.

### Discussion 2.1 - "Are you able to extract the LLP proper decay length from the exponential fit? Does it agree with the expectation?"

Yes. The fit works precisely *because* the distribution is exponential: taking the logarithm turns
$N(t)=e^{-t/\tau}$ into a straight line, so a linear fit on a `Logy` plot is a fit of the lifetime, and the two
samples fall on top of each other as parallel-looking straight lines with different slopes. After the cm $\to$ mm
conversion the extracted values reproduce the sample names: `MC_40_1000` gives $\approx 1000$ mm (1 m) and
`MC_40_10000` gives $\approx 10000$ mm (10 m).

*Pitfall:* dropping the factor 10 makes the answer 10x "wrong" even though the fit itself is perfect, so the
first thing to check when the number looks off is the unit of the fitted variable, not the fit.

---

## Ex2 - "You will implement this equation in the function `weight_calc`"

$$w = \left(\frac{\tau_{old}}{\tau_{new}}\right)^2 \exp\left[(t_1+t_2)\left(\frac{1}{\tau_{old}}-\frac{1}{\tau_{new}}\right)\right]$$

```python
def weight_calc(llp_ct, new_ctau, old_ctau):
    denominator = np.exp(-1.0*llp_ct/old_ctau)/old_ctau**2
    numerator = np.exp(-1.0*llp_ct/new_ctau)/new_ctau**2
    weight = numerator/denominator

    return weight
```

**Explanation.** Instead of coding the formula literally, we code it as *pdf(new) / pdf(old)*, which is
impossible to get wrong algebraically and is exactly how `weight_calc` is written in `lib/helper.py`:

$$p(t_1,t_2\,|\,\tau)=\frac{1}{\tau^2}e^{-t_1/\tau}e^{-t_2/\tau}=\frac{1}{\tau^2}e^{-(t_1+t_2)/\tau}
\;\;\Rightarrow\;\; w=\frac{p(t_1+t_2|\tau_{new})}{p(t_1+t_2|\tau_{old})}.$$

Because the two LLPs decay independently, only their **sum** enters (`llp_ct = t1 + t2`), and the
$(\tau_{old}/\tau_{new})^2$ prefactor is just the ratio of the two $1/\tau^2$ normalisations. `llp_ct`,
`new_ctau`, `old_ctau` must all be in the **same** length unit (cm in this notebook).

### Ex2 (reweight 10 m to 1 m)

```python
        llp_ct = np.sum(gLLP_ctau[k], axis = 1)   # t1 + t2 in cm, already in the skeleton

        # Fill in the arguments in weight_calc()
        # reweight the 10 m sample (old_ctau = 1000 cm) down to 1 m (new_ctau = 100 cm)
        weight = weight_calc(llp_ct, 100.0, 1000.0)

        for j in range(Nevents_todraw): h[k].Fill(gLLP_ctau[k][j][0], weight[j])
        ...
        print("sample: {}, fitted ctau [mm]: {}".format(k, -10./fit.GetParameter(1)))
```

**Explanation.** `np.sum(gLLP_ctau['MC_40_10000'], axis=1)` adds the proper decay lengths of the two LLPs in
each event (axis 1 = the LLP axis, not the event axis). "10 m" and "1 m" must be entered as **1000 cm** and
**100 cm** because the argument order is `(llp_ct, new_ctau, old_ctau)` and `gLLP_ctau` is in cm; `new` comes
before `old`. The weights are then passed as the second argument of `TH1D.Fill`, so the 10 m sample is redrawn
as if it had been generated with a 1 m lifetime. Validation: the reweighted histogram agrees with the genuine
1 m sample and both `expo` fits print ~1000 mm.

### Discussion 2.2 - "Why do we reweight from a larger lifetime to a smaller lifetime? What happens when we do the opposite?"

With $\tau_{new}<\tau_{old}$ the factor $\left(\frac{1}{\tau_{old}}-\frac{1}{\tau_{new}}\right)$ is **negative**,
so $w\propto e^{-t/\tau_{eff}}$: weights only *suppress* late decays and are largest for small $t$, and small-$t$
events are exactly what a long-lived sample is full of. Every event still contributes, the weights are bounded,
and the statistical power of the original sample is preserved.

Reweighting in the opposite direction makes the exponent **positive**: $w\propto e^{+t/\tau_{eff}}$. The input
sample has very few large-$t$ events (they are exponentially rare), so a handful of them receive enormous weights
and dominate the histogram. The mean can look right while the *effective* number of entries collapses
($N_{eff} = (\sum w)^2/\sum w^2 \ll N$), i.e. the result has a huge statistical uncertainty and is unstable under
which events happen to be in the tail. If $\tau_{new}\gg\tau_{old}$ there may be no long-decay events at all, and
the reweighted shape is then pure extrapolation rather than a measurement. This is why the samples generated at
0.1/1/10/100 m are only ever reweighted **down** to the intermediate values 0.5/5/50 m.

---

## Ex3 - "You will first define the geometric decay region of the muon system"

```python
for k, T in tree.items():

     # define endcap muon system region
    sel_ms = np.abs(T['gLLP_eta'].array()) < 2.1
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_z'].array()) < 1050)
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_z'].array()) > 568)
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_r'].array()) < 280)

    print(k, np.count_nonzero(np.sum(sel_ms, axis = 1))/len(sel_ms))
```

**Explanation of the numbers.** The region is the endcap muon system drawn in Fig. 4.1.1 (p. 141) of the Muon
TDR, and it is the same geometry the course uses in `lib/util.py` (`drawCSCz`/`drawCSCsteel`: the grey bands are
the steel absorbers - `preME11`, `ME11_12`, `ME12_2`, `ME2_3`, `ME3_4`, `beyond` - and the CSC chambers sit in the
gaps between them, from the front edge of ME1/1 at |z| = 568 cm to the back of ME4 at |z| = 1050 cm):

| cut | value | meaning |
|---|---|---|
| `abs(gLLP_eta) < 2.1` | - | direction inside the CSC coverage |
| `abs(z) > 568` cm | inner edge | ME1/1, the first endcap station, starts at 568 cm (end of the `preME11` steel) |
| `abs(z) < 1050` cm | outer edge | ME4, the last endcap station, ends at 1050 cm (`beyond` = outside CMS starts there) |
| `abs(r) < 280` cm | radial | size of the endcap disks |

All vertex branches are in cm, so the constants are used directly. `sel_ms` is a boolean array of shape
`(events, 2 LLPs)`; `np.sum(sel_ms, axis=1)` counts how many LLPs of the event decay inside, and
`np.count_nonzero(...)/len(sel_ms)` is the fraction of events with **at least one** such decay = the geometric
acceptance. Expected size: **~15-20 % for the 1 m and 10 m samples**, smaller for 0.1 m (decay happens earlier,
in the calorimeters) and for 100 m (the LLP crosses CMS undecayed).

*Note:* the region is a definition, not a measurement - a looser box (e.g. |eta| < 2.4, 400 < |z| < 1100,
r < 695 cm, which also accepts barrel-like decays) gives the same physics with a somewhat larger probability.
What matters is that the *same* `sel_ms` is used in the acceptance scan below.

---

## Acceptance vs. ctau

```python
    T = tree['MC_40_{}'.format(ct_source)]

    ################## ANSWER #################
    # geometric endcap muon system region, identical to Ex3 (vertices are in cm)
    sel_ms = np.abs(T['gLLP_eta'].array()) < 2.1
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_z'].array()) < 1050)
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_z'].array()) > 568)
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_r'].array()) < 280)

    # ctau_list (ct) and ct_source are in mm, gLLP_ctau is in cm -> divide by 10
    weight = weight_calc(np.sum(gLLP_ctau['MC_40_{}'.format(ct_source)], axis = 1), ct/10.0, ct_source/10.0)
    #################################################

    denominator = np.sum(weight)
    numerator = np.sum(weight[np.sum(sel_ms, axis = 1)>0])
    acceptance.append(numerator/denominator)
```

**Explanation.** For each target lifetime `ct` (in mm) the loop above already chose `ct_source`, the shortest
generated lifetime that is **not smaller** than `ct`, so we always reweight downwards (Discussion 2.2). The
acceptance is the weighted fraction of events with a decay in the muon system: `np.sum(weight)` is the weighted
event count and `np.sum(weight[mask])` is the same for events with at least one in-region LLP. The crucial step is the
**mm $\to$ cm conversion**: `ctau_list` and `OLD_CTAU` are in mm (they mirror the file names) while
`gLLP_ctau` is in cm, hence `ct/10.0` and `ct_source/10.0`. Skipping the `/10.` silently returns the acceptance
of a lifetime 10x shorter - a smooth, plausible-looking, wrong curve.

The final cell plots `ctau_list` (converted to m) vs `acceptance` on a log-log `TGraph`.

### Discussion 2.3 - "Do you understand the shape? Why is the peak at 2 m, much smaller than the distance between the muon detectors and the interaction point?"

The curve must vanish at both ends:

- $c\tau \ll 1$ m: the LLP decays promptly, typically inside the tracker/calorimeters, **before** the muon
  system $\to$ acceptance $\to 0$;
- $c\tau \gg 10$ m: the LLP is effectively stable over the size of CMS and leaves without decaying $\to$
  acceptance $\to 0$;
- in between, the acceptance peaks where the exponential covers the muon-system shell most efficiently.

The peak sits at $c\tau \approx 2$ m rather than at the 6-10 m distance from the IP to the muon system because
the exponential decay happens in the **lab** with the *boosted* length: the detector sees
$\langle L\rangle = \beta\gamma\, c\tau$. These LLPs come from $H\to SS$ with $m_H = 125$ GeV and $m_S = 40$ GeV,
so each LLP carries tens of GeV and $\beta\gamma$ is a few ($\sim 3$-5). Thus $c\tau = 2$ m corresponds to a mean
lab decay length of roughly 6-10 m, exactly the radius of the muon system, and that is where the maximum must
land. Equivalently, the shell at radius $R$ is "illuminated" by an exponential of length $\beta\gamma c\tau$, so
the maximum is at $c\tau \simeq R/\langle\beta\gamma\rangle$, not at $R$. (The finite thickness of the shell and
the spread of $\beta\gamma$ broaden and flatten the peak.) This lifetime dependence is exactly why LLP searches in
different sub-detectors - tracker, calorimeter, muon system - are complementary.
