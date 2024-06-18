## TransductionRHS_v5

\IHC\TransductionRHS_v5.m

https://github.com/evanwporter/cochlea-nerve/blob/19159e7c3610c42b41953590521d990e9c226f77/IHC/TransductionRHS_v5.m#L110-L112

This section determines whether based on membrance potential (Vm) 
channels are open or closed or blocked.
Has an element of randomness (monte carlo)

Calcium current is calculated with GHK

Calcium_concentration uses the calcium current to update the concentration 
of calcium at vesicles (C_vesicle) and channels (C).

Calculate neuro transmitter release using `NTdynamicsRHS_v5_core`

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