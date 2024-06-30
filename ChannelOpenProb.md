# Channel Open Probability

https://github.com/evanwporter/cochlea-nerve/blob/4e8b4f18f20782bfc39d88589db9f4e04dbcf507/IHC/Transduction_v4.m#L390-L396

`z` is the parameter being optimized for

z contains

* Calcium Channel States (# open or blocked)
    -> Affects calcium influx.
* Calcium Concentration (C_vesicles)
* Calcium Current (I)
* And all 4 from `NTdynamicsRHS_v5_core`

# Concentration


$u2 = \left( \frac{r^{2-d}}{4D\pi^{d/2}} \right) \times (\sqrt{\pi} \cdot \text{erfc}(r / \sqrt{4Dt}))$
