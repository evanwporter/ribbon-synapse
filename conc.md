

obj.u2((m+1) * obj.dt) - obj.u2(m * obj.dt)

$u(t) = \frac{r^{2-d}}{4 D \pi^{d/2}}  \sqrt{\pi} \exp{\frac{r}{\sqrt{4 D t}}}$

$E = u(t + dt) - u(t)$

$E$ then is a vector that looks like $\{E_1, E_2, E_3, \ldots, E_n\}$

The code calculates concentration by taking the dot product of the beginning of $E$ and the end of the current array $J$.

To illustrate what I mean, let $J = \{ J_1, J_2, J_3, \ldots, J_n \}$.

It then multiples the two array as follows:

$ J_{n} E_1 + J_{n - 1} E_2 + J_{n - 2} E_3 + \ldots + J_{3} E_{n-2} + J_2 E_{n-1} + J_1 E_{n} $

This is the interesting part, based on what the specified margin of error tolerance is specified the algorithm only takes the last N values of J (and first N values of E).

This is what doesn't make sense to me, why we flip E (I believe he does this because the initial values of E have a larger effect than the later values). Also how he derived U and why he calculated E (change in U over a single time step). The paper doesn't mention this at all. 

Because of these issues I came up with my own method.


$C(t, r) = \int_0^t \frac{I(t')}{(4\pi D (t - t'))^{3/2}} \exp\left(-\frac{|r|^2}{4D (t - t')}\right) \, dt'$

To estimate Concentration I used Rhiemann sums, ie: rectangle method, for solving integrals that I learned about in high school.