# Synthetic coupled dynamical systems

## Three coupled logistic maps

For this example, we'll consider a system of three coupled logistic
maps (``x``, ``y`` and ``z``), where the dynamical influence goes in the directions ``z \rightarrow x`` and ``z \rightarrow y``.

The system is given by the following difference equations:

```math
x(t+1) = [x(t)(r_1 - r_1 x(t) - z(t) + \sigma_{x(t)} \eta_{x}(t))] \mod 1 \\
y(t+1) = [y(t)(r_2 - r_2 y(t) - z(t) + \sigma_{y(t)} \eta_{y}(t))] \mod 1 \\
z(t+1) = [z(t)(r_3 - r_3 z(t) + \sigma_{z(t)} \eta_{z}(t))] \mod 1
```

Dynamical noise may be added to each of the dynamical variables by tuning the parameters `σz`, `σx` and `σz`. Default values for the parameters `r₁`, `r₂` and `r₃` are set such that the system exhibits chaotic behaviour, with `r₁ = r₂ = r₃ = 4`.


## Where has the system been used?
This system was used in [Runge (2018)](https://aip.scitation.org/doi/abs/10.1063/1.5025050), where he
discusses the theoretical assumptions behind causal network reconstructions and practical estimators of such networks.


## Represent as a DiscreteDynamicalSystem

```@setup logistic3
using CausalityTools
using DynamicalSystems
using Plots
using StaticArrays
using Statistics
```

We first define the equations of motion.

```@example logistic3
function eom_logistic3(u, p, t)
    r₁, r₂, r₃, σx, σy, σz = (p...,)
    x, y, z = (u...,)

    # Independent dynamical noise for each variable.
    ηx = rand()
    ηy = rand()
    ηz = rand()

    dx = (x*(r₁ - r₁*x - z + σx*ηx)) % 1
    dy = (y*(r₂ - r₂*y - z + σy*ηy)) % 1
    dz = (z*(r₃ - r₃*z + σz*ηz)) % 1
    return SVector{3}(dx, dy, dz)
end
nothing #hide
```

To make things easier to use, we create function that generates a
DiscreteDynamicalSystem instance for any set of parameters `r₁`, `r₂` and `r₃`, initial condition `u₀`, and dynamical noise levels `σx`, `σy` and `σz`.

```@example logistic3
function logistic3(;u₀ = rand(3), r₁ = 4, r₂ = 4, r₃ = 4,
                    σx = 0.05, σy = 0.05, σz = 0.05)
    p = [r₁, r₂, r₃, σx, σy, σz]
    DiscreteDynamicalSystem(eom_logistic3, u₀, p)
end
nothing #hide
```

An example realization of the system is:

```@example logistic3
s = logistic3()
orbit = trajectory(s, 100)
x, y, z = orbit[:, 1], orbit[:, 2], orbit[:, 3]
plot(x, label = "x", lc = :black)
plot!(y, label = "y", lc = :red)
plot!(z, label = "z", lc = :blue)
xlabel!("Time step"); ylabel!("Value")
savefig("logistic3.svg") #hide
nothing #hide
```

![](logistic3.svg)

## Predefined system
This system is predefined in `CausalityTools.Systems`, and can be initialized using the `logistic3` function.


## References

Runge, Jakob. Causal network reconstruction from time series: From theoretical assumptions to practical estimation, Chaos 28, 075310 (2018); doi: 10.1063/1.5025050. [https://aip.scitation.org/doi/abs/10.1063/1.5025050](https://aip.scitation.org/doi/abs/10.1063/1.5025050)
