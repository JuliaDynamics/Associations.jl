# Synthetic coupled dynamical systems

## Nonlinear response of two periodic drivers
For this example, we'll consider a where the response ``X`` is a highly
nonlinear combination of two drivers ``Y`` and ``Z``.

The system is given by the following difference equations:


```math
\begin{aligned}
x(t+1) &= \dfrac{y(t)(18y(t) - 27y(t)^2 + 10)}{2} + z(t)(1-z(t)) + \eta_x \\
y(t+1) &= \dfrac{1 - \dfrac{\\cos(2 \pi)}{\omega y}t}{2} + \eta_y \\
z(t+1) &= \dfrac{1 - \dfrac{\\sin(2 \pi)}{\omega z}t}{2} + \eta_z,
\end{aligned}
```
where ``η_x``, ``η_y`` and ``η_z`` are gaussian noise terms with mean 0 and standard deviations ``σ_x``, ``σ_y`` and ``σ_z``.


## Where has the system been used?
This system was used in [Verdes (2005)](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.72.026222) where
he studies a nonparametric test for causality in weakly coupled systems.


## Represent as a DiscreteDynamicalSystem

```@setup verdes
using CausalityTools
using DynamicalSystems
using Plots
using StaticArrays
using Statistics
using Distributions
```

We first define the equations of motion.

```@example verdes
function eom_verdes(u, p, t)
    x, y, z = (u...,)
    ωy, ωz, σx, σy, σz = (p...,)

    ηx = σx == 0 ? 0 : rand(Normal(0, σx))
    ηy = σy == 0 ? 0 : rand(Normal(0, σy))
    ηz = σz == 0 ? 0 : rand(Normal(0, σz))

    dx = y*(18y - 27y^2 + 10)/2 + z*(1-z) + ηx
    dy = (1 - cos((2*pi/ωy) * t))/2 + ηy
    dz = (1 - sin((2*pi/ωz) * t))/2 + ηz
    return SVector{3}(dx, dy, dz)
end
nothing #hide
```

To make things easier to use, we create function that generates a
DiscreteDynamicalSystem instance for any set of parameters `r₁`, `r₂` and `r₃`, initial condition `u₀`, and dynamical noise levels `σx`, `σy` and `σz`.

```@example verdes
function verdes(;u₀ = rand(3), ωy = 315, ωz = 80,
                σx = 0.01, σy = 0.01, σz = 0.01)
    p = [ωy, ωz, σx, σy, σz]
    DiscreteDynamicalSystem(eom_verdes, u₀, p)
end
nothing #hide
```

An example realization of the system is:

```@example verdes
s = verdes()
orbit = trajectory(s, 500)
x, y, z = orbit[:, 1], orbit[:, 2], orbit[:, 3]
plot(x, label = "x", lc = :black)
plot!(y, label = "y", lc = :red)
plot!(z, label = "z", lc = :blue)
xlabel!("Time step"); ylabel!("Value")
savefig("verdes.svg") #hide
nothing #hide
```

![](verdes.svg)

## Predefined system
This system is predefined in `CausalityTools.Systems`, and can be initialized using the `verdes` function.


## References
Verdes, P. F. "Assessing causality from multivariate time series." Physical Review E 72.2 (2005): 026222. [https://journals.aps.org/pre/abstract/10.1103/PhysRevE.72.026222](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.72.026222)
