# Synthetic coupled dynamical systems

## Two unidirectionally coupled logistic maps

For this example, we'll consider a unidirectionally coupled system consisting
of two logistic maps, given by the difference equations

```math
\begin{aligned}
dx &= r_1 x(1 - x) \\
dy &= r_2 f(x,y)(1 - f(x,y)),
\end{aligned}
```

with

```math
\begin{aligned}
f(x,y) = \dfrac{y + \frac{c(x + \sigma \xi )}{2}}{1 + \frac{c}{2}(1+ \sigma )}
\end{aligned}
```

The parameter `c` controls how strong the dynamical forcing is. If `σ > 0`,
dynamical noise masking the influence of  ``x`` on ``y``, equivalent to
``\sigma \cdot \\xi``, is added at each iteration. Here, ``\xi`` is a draw from
a flat distribution on ``[0, 1]``. Thus, setting `σ = 0.05` is equivalent to
add dynamical noise corresponding to a maximum of ``5 \%`` of the possible
range of values of the logistic map.

## Where has the system been used?
This system was used in [Diego et al. (2018)](https://ui.adsabs.harvard.edu/#abs/2018arXiv181101677D)
to study the performance of the transfer operator grid transfer
entropy estimator.

## Represent as a DiscreteDynamicalSystem

```@setup logistic2
using CausalityTools
using DynamicalSystems
using Plots
using Statistics
```

We first define the equations of motion.

```@example logistic2
function eom_logistic2(dx, x, p, n)
    c, r₁, r₂, σ = (p...,)
    ξ = rand() # random number from flat distribution on [0, 1]
    x, y = x[1], x[2]
    f_xy = (y +  (c*(x + σ*ξ)/2) ) / (1 + (c/2)*(1+σ))

    dx[1] = r₁ * x * (1 - x)
    dx[2] = r₂ * (f_xy) * (1 - f_xy)
    return
end
nothing #hide
```

To make things easier to use, we create function that generates a
DiscreteDynamicalSystem instance for any set of parameters `r₁` and `r₂`,
coupling strength `c`, initial condition `u₀` and dynamical noise level `σ`.
Selecting parameter values on `[3.6, 4.0]` yield mostly chaotic realizations of
the maps, so we set the default to some random values on this interval.

```@example logistic2
function logistic2(;u₀ = rand(2), c = 0.0, r₁ = 3.66, r₂ = 3.77, σ = 0.05)
    p = [c, r₁, r₂, σ]
    DiscreteDynamicalSystem(eom_logistic2, u₀, p)
end
nothing #hide
```

By tuning the coupling strength `c`, we may control the strength of the influence
``x`` has on ``y``. Depending on the particular values of `r₁` and `r₂`, the
subsystems become synchronized at different values of `c`. Choosing `c ∈ [0, 2]`
usually still gives some independence between the subsystems.

An example realization of the system when there is no coupling is:

```@example logistic2
s = logistic2(c = 0.0)
orbit = trajectory(s, 100)
x, y = orbit[:, 1], orbit[:, 2]
plot(x, label = "x", lc = :black)
plot!(y, label = "y", lc = :red)
xlabel!("Time step"); ylabel!("Value")
savefig("logistic2_c0.svg") #hide
nothing #hide
```

![](logistic2_c0.svg)

## Predefined system
This system is predefined in `CausalityTools.Systems`, and can be initialized using the `logistic2` function.


## References
Diego, D., Agasøster Haaga, K., & Hannisdal, B. (2018, November 1).
Transfer entropy computation using the Perron-Frobenius operator.
Eprint ArXiv:1811.01677. Retrieved from
[https://ui.adsabs.harvard.edu/#abs/2018arXiv181101677D](https://ui.adsabs.harvard.edu/#abs/2018arXiv181101677D)
