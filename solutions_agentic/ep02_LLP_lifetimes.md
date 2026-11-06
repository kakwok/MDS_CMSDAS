# Episode 02: LLP lifetimes (solutions, agentic)

## Ex0 - LLP proper decay length (cell 7)

```python
gLLP_ctau = {}
for k, T in tree.items():
    gLLP_beta = T['gLLP_beta'].array()
    gLLP_decay_vertex = np.sqrt(T['gLLP_decay_vertex_x'].array()**2 + T['gLLP_decay_vertex_y'].array()**2 + T['gLLP_decay_vertex_z'].array()**2)
    gLLP_gamma = 1.0/np.sqrt(1.0 - gLLP_beta**2)
    gLLP_ctau[k] = gLLP_decay_vertex/(gLLP_beta*gLLP_gamma)
```

d = |decay vertex| is the lab-frame decay length (cm), since the LLP is produced at the origin. Then d = beta*gamma*ctau gives ctau = d/(beta*gamma), with gamma = 1/sqrt(1-beta^2). Units stay in cm, the unit of the branches.

## Ex1 - plot and exponential fit (cell 10)

```python
%%time

Nevents_todraw = 10000
c = rt.TCanvas('c','c', 800, 600)
h = {}

for i, k in enumerate([ 'MC_40_1000', 'MC_40_10000']):
    h[k] = rt.TH1D('', '', 100, 0, 10000)
    h[k].SetXTitle('LLP c#tau [cm]')
    h[k].SetYTitle('fraction of events')
    for j in range(Nevents_todraw): h[k].Fill(gLLP_ctau[k][j][0])
    h[k].Scale(1./h[k].Integral())  #normalize histogram to 1
    h[k].SetLineColor(i+1)
    h[k].Draw('hist same')

    ####### perform exponential fit ######
    h[k].Fit('expo')
    fit = h[k].GetFunction('expo')
    # ROOT 'expo' is exp(p0 + p1*x), so p1 = -1/ctau with ctau in cm
    fitted_slope = fit.GetParameter(1)
    print("sample: {}, fitted ctau [mm]: {}".format(k, -10.0/fitted_slope))
    #######################################

c.SetLogy()
c.Draw()
```

**Explanation (Ex1).** The histogram of gLLP_ctau is normalised to unit area and fitted with ROOT's `expo`, which is exp(p0 + p1*x). On a log-y plot this is a straight line, which checks the exponential shape. The fitted slope p1 equals -1/ctau, in the histogram unit (cm). The printed value is converted to mm with a factor of 10, because the file names are in mm. Pitfall: the 1 m sample has ctau of about 100 cm, which is about the same as the 100 cm bin width here, so the fit is coarse; the 10 m sample is better sampled.

**Discussion 2.1 (answer).** Question, verbatim: "Are you able to extract the LLP proper decay length from the exponential fit? Does it agree with the expectation?"

Yes. The fitted ctau is -10/p1 in mm. For the 1 m (MC_40_1000) and 10 m (MC_40_10000) samples the printed values should be close to 1000 and 10000 mm. The 10 m sample should agree within a few percent. The 1 m sample can deviate more because of the coarse binning and the fit range, so the agreement is at the level of the bin width and the fit uncertainty. I did not execute the notebook, so the exact printed values are not quoted here; check them against the expected values.

## Ex2 - weight_calc (cell 13)

```python
'''
arguements:
llp_ctau is the list of sum of proper decay length, which is t1+t2 in the equation
new_ctau is tau_new in the equation
old_ctau is tau_old in the equation
All units in cm

return event weight

'''


def weight_calc(llp_ct, new_ctau, old_ctau):
    denominator = np.exp(-1.0*llp_ct/old_ctau)/old_ctau**2
    numerator = np.exp(-1.0*llp_ct/new_ctau)/new_ctau**2
    weight = numerator/denominator

    return weight
```

**Explanation (Ex2).** The event probability for two independent LLPs is p(t1,t2|tau) = (1/tau^2) exp(-(t1+t2)/tau). The weight from a sample with tau_old to tau_new is the ratio p(new)/p(old). Dividing by the old density and multiplying by the new density gives w = (tau_old/tau_new)^2 exp[(t1+t2)(1/tau_old - 1/tau_new)], which is the formula in the lesson. The code computes the same thing as numerator/denominator. Units: llp_ct, new_ctau and old_ctau must all be in the same unit, which is cm here. The lib/helper.py version does the same job, but this notebook keeps the function inline.

## Ex2 - reweight 10 m to 1 m and validate (cell 16)

```python
%%time
c = rt.TCanvas('c','c', 800, 600)
h = {}
Nevents_todraw = 10000

for i, k in enumerate(['MC_40_1000', 'MC_40_10000']):
    h[k] = rt.TH1D('', '', 100, 0, 1000)
    h[k].SetXTitle('LLP c#tau [cm]')
    h[k].SetYTitle('fraction of events')

    if k == 'MC_40_1000':
        for j in range(Nevents_todraw): h[k].Fill(gLLP_ctau[k][j][0])
    else:
        llp_ct = np.sum(gLLP_ctau[k], axis = 1)

        # 10 m (1000 cm) -> 1 m (100 cm), both in cm
        weight = weight_calc(llp_ct, 100.0, 1000.0)

        for j in range(Nevents_todraw): h[k].Fill(gLLP_ctau[k][j][0], weight[j])

    h[k].Scale(1./h[k].Integral())
    h[k].SetLineColor(i+1)
    h[k].Draw('hist same')


    ####### perform exponential fit ######
    h[k].Fit('expo')
    fit = h[k].GetFunction('expo')
    # ROOT 'expo' is exp(p0 + p1*x), so p1 = -1/ctau with ctau in cm
    fitted_slope = fit.GetParameter(1)
    print("sample: {}, fitted ctau [mm]: {}".format(k, -10.0/fitted_slope))

    #######################################

c.SetLogy()
c.Draw()
```

**Explanation (Ex2 validation).** The 10 m sample is reweighted to 1 m with weight_calc(llp_ct, 100 cm, 1000 cm). llp_ct is the sum of the two LLP proper decay lengths in each event, in cm. The reweighted 10 m histogram should match the directly generated 1 m histogram, and its fitted ctau should be about 1000 mm. The fit is printed in mm, using the same factor of 10 as Ex1.

**Discussion 2.2 (answer).** Question, verbatim: "Why do we reweight from a larger lifetime to a smaller lifetime? What happens when we do the opposite (reweight from smaller lifetime to larger lifetime)?"

Reweighting from a larger to a smaller lifetime is safe. The old sample has events everywhere the new distribution needs them, because the small-lifetime region is well populated in the long-lived sample. The weight exp[(t1+t2)(1/tau_old - 1/tau_new)] is then bounded and close to 1 for the events that matter. Going the other way, the new distribution is concentrated at large t, where the old, short-lifetime sample has almost no events. Its weights become very large, a few events dominate, the effective sample size drops, and the statistical uncertainty grows. The ctau_list code in Ex3 follows this rule: it picks the smallest generated sample with tau_old >= tau_new, and the extrapolation beyond 100 m is the one place where the rule is broken.

## Ex3 - endcap muon decay region (cell 19)

```python
%%time

for k, T in tree.items():

     # define endcap muon system region
    sel_ms = np.abs(T['gLLP_eta'].array()) < 2.4
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_z'].array()) < 1075)
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_z'].array()) > 568)
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_r'].array()) < 700)

    print(k, np.count_nonzero(np.sum(sel_ms, axis = 1))/len(sel_ms))
```

**Explanation (Ex3 region).** The endcap muon region is set from the CSC layout in lib/util.py (drawRZ) and the CSC z-plot used in MDS_reconstruction. Cuts (cm, cm, cm, and eta):
- |eta| < 2.4: muon system acceptance.
- 568 < |z| < 1075: from ME1/1 (z about 568 to 632 cm) to the last station, with the 1075 cm limit as the outer edge.
- |r| < 700 cm: outer radius of the endcap chambers (ME4/2 and ME1/3 reach about 700 cm).

The printed value is the fraction of events with at least one LLP decaying in this volume. The target is about 15 to 20 percent for the 1 m and 10 m samples. I did not run the notebook, so the value is not confirmed. Pitfall: the z lower bound must be an absolute value, because the endcaps are at both +z and -z.

## Ex3 - acceptance vs ctau (cell 22)

```python
%%time

ctau_list  = [100, 300, 500, 1000, 2000, 5000, 10000, 50000, 100000, 500000] #mm

acceptance = []
for ct in ctau_list:
    # find which sample to use for reweighting
    ct_source = 0
    for i in OLD_CTAU:
        if ct <= i:
            ct_source = i
            break
    if ct_source == 0: ct_source = OLD_CTAU[-1]
    print(ct, ct_source) #verify the correct sample is chosen

    T = tree['MC_40_{}'.format(ct_source)]

    ################## FILL IN HERE #################
    sel_ms = np.abs(T['gLLP_eta'].array()) < 2.4
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_z'].array()) < 1075)
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_z'].array()) > 568)
    sel_ms = np.logical_and(sel_ms, np.abs(T['gLLP_decay_vertex_r'].array()) < 700)
    # sum of the two LLP proper decay lengths in cm, from the source sample
    llp_ct = np.sum(gLLP_ctau['MC_40_{}'.format(ct_source)], axis=1)
    # reweight from the source lifetime (mm -> cm) to the target lifetime (mm -> cm)
    weight = weight_calc(llp_ct, ct/10.0, ct_source/10.0)

    #################################################

    denominator = np.sum(weight)
    numerator = np.sum(weight[np.sum(sel_ms, axis = 1)>0])
    acceptance.append(numerator/denominator)

print(ctau_list, acceptance)
```

**Explanation (Ex3 acceptance).** For each target ctau the code picks the smallest generated sample with ctau_source >= ctau, reweights its events with weight_calc, and computes acceptance = (sum of weights with at least one LLP in the endcap region) / (sum of all weights). The units are converted once: the ntuple quantities are in cm, the ctau_list is in mm, so the target and source are divided by 10 before the call. The region cut is the same one as in the cell above, so the two cells cannot disagree.

**Q: Do you understand the shape? (answer).** The notebook question is "Q: Do you understand the shape?".

Yes. Acceptance is small at short ctau because most LLPs decay before reaching the muon system, which is at 4 to 10 m from the IP. Acceptance falls at long ctau because most LLPs fly past the detector before they decay. The maximum sits where the lab-frame decay length, beta*gamma*ctau, is comparable to the distance to the endcap.

**Discussion 2.3 (answer).** Question, verbatim: "Do you understand the shape of the geometric acceptance vs. LLP mean proper decay lengths plot? Why is the peak at 2m, much smaller than the distance between the muon detectors and the interaction point?"

The peak is at ctau of about 2 m, not at the detector distance, because the quantity that must match the distance is the lab-frame decay length beta*gamma*ctau, not ctau. The LLPs in this sample are boosted (beta*gamma of a few), so a ctau of about 2 m gives lab decay lengths of several meters, which is the detector distance scale. Acceptance is a product of two factors: the chance of surviving to the endcap, which is exp(-d/lambda), and the chance of decaying inside the endcap volume, which is 1 - exp(-Delta/lambda). Their product peaks at lambda of order d and Delta, not at the full distance d. The exact value depends on the boost distribution in the sample, which can be checked with gLLP_beta.

