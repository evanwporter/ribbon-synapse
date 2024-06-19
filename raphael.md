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

## NTdynamics

Calculates neurotransmitter change in states over time (dt).

Returns the following:
https://github.com/evanwporter/cochlea-nerve/blob/cc5499247e8e09865531b8ef920a54d3ea39d4cc/IHC/NTdynamicsRHS_v5_core.m#L41-L45


## TransductionRHS_v5

https://github.com/evanwporter/cochlea-nerve/blob/19159e7c3610c42b41953590521d990e9c226f77/IHC/TransductionRHS_v5.m#L1

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

## Transduction_v4

https://github.com/evanwporter/cochlea-nerve/blob/a9d084b5d8723044151d401e4a037e517410a432/IHC/Transduction_v4.m#L1-L2

There's a lot going on here but this is where he does the ODE stuff. He uses a function called odeEuler to solve the ODE.

https://github.com/evanwporter/cochlea-nerve/blob/3146900feea2064af26d5eec043a1f72d12819ad/Wrapper/Miscellaneous/odeEuler.m#L1-L2

odeEuler is an extrememly simple ODE solver. Solves the following equation

$y_{n+1} = y_n + \delta t * f(t_n, y_n)$

But basically he uses both `TransductionRHS_v5` and `NTdynamicsRHS_v5` as the ODEFUNC. This is probably the file to look at.

He also used a custom function odeWrapper, which is just an interface for handling the results.

https://github.com/evanwporter/cochlea-nerve/blob/72247ca25071fc91b1c5b95c4a53867c2a4309cf/IHC/Transduction_v4.m#L382-L387

It worth noting that if we plan to use `odeWrapper` in any significant way, then we may need to verify how it works since the author express some uncertainty over the mathematical components:

https://github.com/evanwporter/cochlea-nerve/blob/72247ca25071fc91b1c5b95c4a53867c2a4309cf/externals/ode-wrapper/odeWrapper.m#L1-L4





<!-- |-> NTdynamicsRHS_v5_core
|-> TransductionRHS_v5
|-> Transduction_v4
|-> Transduction_v4_multi
|-> Synapse
|-> setSynapseResult
|-> handleResult
|-> handleSynapseResult
|-> ANT / ANT_clamp
|-> oneANTjob / runANTmulti
|-> ex4_2_synapse -->