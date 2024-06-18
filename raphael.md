This is a massive program with over 2000 files! The part I've looked at make up only a small fraction of this program. So I couldn't tell you what a lot of the code is doing. The parts I've focused and the things I could tell you about are the parts relating to the research done in the Raphael Lab. More specifically the modeling done to explain the interactions between the Ribbon Synapse, the calcium curren and the neurotransmitter release.

One of my takeaways from Section 2.9 in the paper was that the release of the neurotransmitter is a function of Ca^{2+}.

We can start at the file [`RibbonSynapse_v4.m`](https://github.com/evanwporter/cochlea-nerve/blob/769ac0d265121edeef3496b53732817e2f28ef66/IHC/RibbonSynapse_v4.m#L66).

## Ribbon Synapse

This class models the structure of the ribbon synapse. This class does two things:

1) **Models calcium channels and vesicles**: Randomly distributes channel and vesicle according to some shape (ie: `uniform ring`)

This is important for modeling and visualizing. See below

2) **Calculate Rate of Release of Neurotransmitters**: Uses Hill-Langmuir equation based on calcium concentration

https://github.com/evanwporter/cochlea-nerve/blob/2a05a7d3f15e9140e37227bb9964407e91b70c09/IHC/RibbonSynapse_v4.m#L98-L107

TransmitterRelease is called by `NTdynamicsRHS_v5_core.m`.

# NTdynamics

Calculates neurotransmitter change in states over time (dt).

Returns the following:
https://github.com/evanwporter/cochlea-nerve/blob/cc5499247e8e09865531b8ef920a54d3ea39d4cc/IHC/NTdynamicsRHS_v5_core.m#L41-L45


## TransductionRHS_v5

https://github.com/evanwporter/cochlea-nerve/blob/19159e7c3610c42b41953590521d990e9c226f77/IHC/TransductionRHS_v5.m#L110-L112

First part of TransductionRHS determines whether based on membrance potential (`Vm`) channels are open or closed or blocked. It has an element of randomness (monte carlo)

Second part calculates calcium current using GHK. Then it calculates the calcium concentration useing the calcium current of calcium at vesicles (C_vesicle) and channels (C).

Calculate neuro transmitter release using `NTdynamicsRHS_v5_core`

Finally combines all this information to find the change over time in the following state variables

* Calcium Channel States (# open or blocked)
    -> Affects calcium influx.
* Calcium Concentration (C_vesicles)
* Calcium Current (I)
* And all 4 from `NTdynamicsRHS_v5_core`

Combines these into `dz`.

https://github.com/evanwporter/cochlea-nerve/blob/19159e7c3610c42b41953590521d990e9c226f77/IHC/NTdynamicsRHS_v5_core.m#L37


|-> NTdynamicsRHS_v5_core
|-> TransductionRHS_v5
|-> Transduction_v4
|-> Transduction_v4_multi
|-> Synapse
|-> setSynapseResult
|-> handleResult
|-> handleSynapseResult
|-> ANT / ANT_clamp
|-> oneANTjob / runANTmulti
|-> ex4_2_synapse