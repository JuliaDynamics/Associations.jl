# Synthetic coupled dynamical systems

## Bivariate AR1 system
For this example, we'll consider a bivariate, order one autoregressive model
consisting of variable ``x`` and ``y``, where ``x \rightarrow y``, given
by the difference equations:

```math
\begin{aligned}
x(t+1) &= a_1 x(t) + \xi_1(t) \\
y(t+1) &= b_1 y(t) + c_1 x(t)+ \xi_2(t)
\end{aligned}
```

where the parameter ``c_1`` controls how strong the dynamical forcing from
``x`` to ``y`` is, and
``\xi_1`` and ``\xi_2`` are dynamical noise terms with zero mean and standard
deviations ``\sigma``.

## Where has the system been used?
This system was investigated by [Paluš et al. (2018)](http://doi.org/10.1063/1.5019944).


## Represent as a DiscreteDynamicalSystem

```@setup ar1
using CausalityTools
using DynamicalSystems
using Plots
using Statistics
using StaticArrays
using Distributions
```

We first define the equations of motion.

```@example ar1
function eom_ar1(x, p, n)
    a₁, b₁, c₁, σ = (p...,)
    x, y = (x...,)
    ξ₁ = rand(Normal(0, σ))
    ξ₂ = rand(Normal(0, σ))

    dx = a₁*x + ξ₁
    dy = b₁*y + c₁*x + ξ₂
    return SVector{2}(dx, dy)
end
nothing #hide
```

To make things easier to use, we create function that generates a
DiscreteDynamicalSystem instance for any
coupling strength `c` and initial condition `u₀`.


```@example ar1
function ar1(;uᵢ = rand(2), a₁ = 0.90693, b₁ = 0.40693, c₁ = 0.5, σ = 0.40662)
    p = [a₁, b₁, c₁, σ]
    return DiscreteDynamicalSystem(eom_ar1, uᵢ, p)
end
nothing #hide
```

By tuning the coupling strength `c₁`, we may control the strength of the influence ``x`` has on ``y``.
An example realization of the system when the coupling strength is `c₁ = 0.5` is:

```@example ar1
s = ar1(c₁ = 0.5)
orbit = trajectory(s, 100)
x, y = orbit[:, 1], orbit[:, 2]
plot(x, label = "x", lc = :black)
plot!(y, label = "y", lc = :red)
xlabel!("Time step"); ylabel!("Value")
savefig("ar1.svg") #hide
nothing #hide
```

![](ar1.svg)

## Predefined system
This system is predefined in `CausalityTools.Systems`, and can be initialized using the `ar1` function.

## References
Paluš, M., Krakovská, A., Jakubík, J., & Chvosteková, M. (2018). Causality,
dynamical systems and the arrow of time. Chaos: An Interdisciplinary Journal of
Nonlinear Science, 28(7), 075307. [http://doi.org/10.1063/1.5019944](http://doi.org/10.1063/1.5019944)
