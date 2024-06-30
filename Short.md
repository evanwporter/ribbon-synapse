# Short Explanation

## What It Does

* Calculates Neurotransmitter release given a voltage

![rrate](https://github.com/evanwporter/ribbon-synapse/assets/115374841/e05fda34-cb92-4942-9b8b-3fcf107b55a8)

* Models Ribbon Synapse

![image](https://github.com/evanwporter/ribbon-synapse/assets/115374841/2b8a071d-d19e-449d-a228-a7f7d00c1db4)

## What it can do

* Calculate Neurotransmitter contents in the vesicle
* Calculate RoC of Neurotransmitter concentration in the synaptic cleft
* Calculate Neurotransmitter Reprocessing rate
* Calculate RoC of proton concentration in the synaptic cleft
* Model
  * RoC of fraction of open channels
  * RoC of fraction of blocked channels
* RoC of calcium current
* RoC of calcium concentration at the vesicles

These things combine into a state vector which it solves as an ode

## How it does it 

* Uses point source diffusion equation to calculate concentration

![image](https://github.com/evanwporter/ribbon-synapse/assets/115374841/2434f9ce-7f23-4597-b4ec-69ec6658627d)

* GHK to calculate current
* Calculate whether channel open or close (element of randomness here)
* Homemade ODE solver that essentially loops through this $y_{n+1} = y_n + \Delta t * f(t_n, y_n)$
* The author chose to model the ribbon synapse a certain way, so I did what he did, but the possibilities for configurations of the rs are unlimited

## Room For Improvement

* `ODE15s` for stiff ODE. Current ODE solver may have issues with accuracy and computatationally be slow since it requires extremely small time steps to be effective.
* Better method of calculating concentration of diffusion
