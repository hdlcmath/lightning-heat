# lightning-heat

These are some basic, self-contained scripts that replicates the main ideas from the paper [The Lightning Method for the Heat Equation](https://link.springer.com/article/10.1007/s10915-026-03222-x).

mhe.m solves a single modified Helmholtz problem, heat.m solves the appropriate number in the Talbot contour and then computes the heat solution, and cumulative_flux.m solves it for various times and computes by trapezoidal rule the cumulative flux.

(Note that in cumulative_flux.m, we actually solve for $\frac{\hat{u}}{s}$, which is computationally more stable and yields the cumulative flux directly as opposed to the flux at a given time, which can be recovered by differentiation.)
