# Concentration Equation (Diffusion)

For each Channel, concentration is calculated at two points of interest:

1. At the vesicle
2. At the channels

It is a function of distance from source channel to point of interest and time.

The equation we are trying to solve here is the single point source diffusion equation:

$\frac{\partial C}{\partial t} = D \nabla^2 C(t, |r|) + I(t) \delta(|r|)$

Where $I$ Is the current through the channel. It is calculated using the GHK equation:

$\Phi = P_{Ca} z^2 \frac{VF^2}{RT} \frac{C_{in} - C_{out} \exp(-zVF/RT)}{1 - exp(-zVF / RT)}$

The author calculated $C(|r|, t)$ as follows

$u(t) = \frac{|r|^{2-d}}{4 D \pi^{d/2}}  \sqrt{\pi} \exp{\frac{|r|}{\sqrt{4 D t}}}$

Where:

* $|r|$ : Distance from point source to point interest
* $D$ : Diffusion constant
* $d$ : Dimensionality (is that a word?). Can be 1, 2, or 3 but in our code this is 3.
* $t$ : 1 x n array of time steps

$E = u(t + dt) - u(t)$

$E$ then is a vector that looks like $\{E_1, E_2, E_3, \ldots, E_n\}$

The code calculates concentration by taking the dot product of the beginning of $E$ and the end of the current array $J$.

To illustrate what I mean, let $J = \{ J_1, J_2, J_3, \ldots, J_n \}$.

It then takes the dot product of both arrays as follows:

$J_{n} E_1 + J_{n - 1} E_2 + J_{n - 2} E_3 + \ldots + J_{3} E_{n-2} + J_2 E_{n-1} + J_1 E_{n}$

To get Concentration $C(|r|, t)$

As a side note since we measure it in three dimensions then the algo multiplies it by 2.

This is the interesting part, based on what the specified margin of error tolerance is specified the algorithm only takes the last N values of J (and first N values of E).

This is what doesn't make sense to me, why we flip E (I believe he does this because the initial values of E have a larger effect than the later values). Also how he derived U and why he calculated E (change in U over a single time step). The paper doesn't mention this at all. 

When N = 1 this is equivalent to

$2 * \frac{J_i}{4 \pi D r} \exp(\frac{r}{\sqrt{4 D t}})$

## Evan's Method
Because of these issues I came up with my own method.

$C(t, r) = \int_0^t \frac{I(t')}{(4\pi D (t - t'))^{3/2}} \exp(-\frac{|r|^2}{4D (t - t')}) \, dt'$

To estimate Concentration I used Rhiemann sums, ie: the rectangle method, for solving integrals that I learned about in high school.

I had to do a lot of online research and thinking to obtain this function. 


To obtain this I used Green's Function

$G(r, t) = \frac{1}{(4 \pi D t)^{3/2}} \exp(-\frac{r^2}{4 D t})$

Next we need to use convolution [^1][^2]:

$(f * g)(t)=\int^{\infty}_{-\infty} f(\tau)g(t - \tau)\$


### Helpful Sources (in no particular order) (yes I know wikipedia is not the best source but it was still helpful so I included it)
## Concentration Comparison Images

Generated using sinusoidal current.

![conc comp](https://github.com/evanwporter/ribbon-synapse/assets/115374841/e196df94-e784-4543-9c71-15a472e93b08)

![conc comp](https://github.com/evanwporter/ribbon-synapse/assets/115374841/983ce530-4a77-475a-95a8-1897176f74db)

![conc comp](https://github.com/evanwporter/ribbon-synapse/assets/115374841/fe1c5f38-91f1-44bb-831f-1b9ca4a03bf9)

![conc comp](https://github.com/evanwporter/ribbon-synapse/assets/115374841/93f772df-0aac-4e6e-a4ea-2e650a5b350c)

![conc comp](https://github.com/evanwporter/ribbon-synapse/assets/115374841/a05ba945-d650-488f-867b-e6a11e87fe78)


## Sources I Found Helpful (in no particular order) (yes I know wikipedia is not the best source but it was still helpful so I included it)

* https://www.physics.uci.edu/~silverma/bseqn/bs/node5.html

* https://en.wikipedia.org/wiki/Green%27s_function

* https://en.wikipedia.org/wiki/Heat_equation

* https://web.pdx.edu/~daescu/mth428_528/Green_functions.pdf

* https://www.physics.uci.edu/~silverma/bseqn/bs/node7.html

* https://en.wikipedia.org/wiki/Fick%27s_laws_of_diffusion

* https://www.rose-hulman.edu/~bryan/lottamath/heatkern.pdf

* https://math.nyu.edu/~tabak/PDEs/The_Heat_Equation.pdf

* https://en.wikipedia.org/wiki/Convection%E2%80%93diffusion_equation

* https://math.stackexchange.com/questions/2588931/green-s-function-for-the-heat-equation

* https://www.damtp.cam.ac.uk/user/dbs26/1BMethods/GreensPDE.pdf


[^1]: https://betterexplained.com/articles/intuitive-convolution/

[^2]: https://developer.nvidia.com/discover/convolution


