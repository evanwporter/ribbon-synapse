# Concentration Equation (Diffusion)

$u(t) = \frac{r^{2-d}}{4 D \pi^{d/2}}  \sqrt{\pi} \exp{\frac{r}{\sqrt{4 D t}}}$

$E = u(t + dt) - u(t)$

$E$ then is a vector that looks like $\{E_1, E_2, E_3, \ldots, E_n\}$

The code calculates concentration by taking the dot product of the beginning of $E$ and the end of the current array $J$.

To illustrate what I mean, let $J = \{ J_1, J_2, J_3, \ldots, J_n \}$.

It then takes the dot product of both arrays as follows:

$ J_{n} E_1 + J_{n - 1} E_2 + J_{n - 2} E_3 + \ldots + J_{3} E_{n-2} + J_2 E_{n-1} + J_1 E_{n} $

As a side note since we measure it in three dimensions then the algo multiplies it by 2.

This is the interesting part, based on what the specified margin of error tolerance is specified the algorithm only takes the last N values of J (and first N values of E).

This is what doesn't make sense to me, why we flip E (I believe he does this because the initial values of E have a larger effect than the later values). Also how he derived U and why he calculated E (change in U over a single time step). The paper doesn't mention this at all. 

When N = 1 this is equivalent to

$2 * \frac{J_i}{4 \pi D r} \exp(\frac{r}{\sqrt{4 D t}})$

## Evan's Method
Because of these issues I came up with my own method.

$C(t, r) = \int_0^t \frac{I(t')}{(4\pi D (t - t'))^{3/2}} \exp\left(-\frac{|r|^2}{4D (t - t')}\right) \, dt'$

To estimate Concentration I used Rhiemann sums, ie: the rectangle method, for solving integrals that I learned about in high school.

I had to do a lot of online research and thinking to obtain this function. 

### Helpful Sources (in no particular order) (yes I know wikipedia is not the best source but it was still helpful)

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

