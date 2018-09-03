# Examples of coupled dynamical systems

The `DynamicalSystems` package implements many famous example of dynamical
systems that have been studied intensively in the literature. However, not many of these systems are coupled in a manner that is well-suited for causality detection algorithms. Such algorithms intend to measure the information flow,
or dynamical coupling strength, between time series. This requires dynamical
systems that are coupled such that we can control the relative dynamical influence between variables.

Numerous such coupled systems appear in the causality detection literature.
Below follows on overview of the different systems, examples of their
usage and literature references.

```@setup s
using CausalityTools
using Plots
using DynamicalSystems
```

# Discrete maps
## Two coupled logistic maps
```@docs
CausalityTools.Systems.logistic2(;u₀ = rand(2), c = 2.0, r₁ = 3.78, r₂ = 3.66)
```

```@example s
# Initialise an instance of the system where `x`
# drives `y` with coupling strength `c = 1.5`.
ds = CausalityTools.Systems.logistic2(c = 1.5)
orbit = trajectory(ds, 100)

# Plot the orbit
t = 1:size(orbit, 1)
plot(t, orbit[:, 1], label = "x")
plot!(t, orbit[:, 2], label = "y")
title!("Two coupled logistic maps, forcing X -> Y")
```

### References
D Diego, KA Haaga, B Hannisdal, in prep. Transfer Entropy computation by Perron-Frobenius operator approximation.

## Three coupled logistic maps
```@docs
CausalityTools.Systems.logistic3(;u₀ = rand(3), r = 4, σx = 0.05, σy = 0.05, σz = 0.05)
```

```@example s
# Create an instance of the discrete dynamical system and
# produce an orbit of 100 points.
ds = CausalityTools.Systems.logistic3()
orbit = trajectory(ds, 100)

# Plot the orbit
t = 1:size(orbit, 1)
plot(t, orbit[:, 1], label = "X")
plot!(t, orbit[:, 2], label = "Y")
plot!(t, orbit[:, 3], label = "Z")
xlabel!("Time step")
ylabel!("Value")
title!("Three coupled logistic maps, Z -> X and Z -> Y")
```

## Four coupled logistic maps

```@docs
CausalityTools.Systems.logistic4(;u₀ = rand(4), r₁ = 3.9, r₂ = 3.6, r₃ = 3.6, r₄ = 3.8, c₂ = 0.4, c₃ = 0.4, c₄ = 0.35)
```

```@example s
# Create an instance of the discrete dynamical system and
# produce an orbit of 100 points.
ds = CausalityTools.Systems.logistic4()
orbit = trajectory(ds, 100)

# Plot the orbit
t = 1:size(orbit, 1)
plot(t, orbit[:, 1], label = "X1")
plot!(t, orbit[:, 2], label = "X2")
plot!(t, orbit[:, 3], label = "X3")
plot!(t, orbit[:, 4], label = "X4")
xlabel!("Time step")
ylabel!("Value")
title!("4 coupled logistic maps, X1 -> X2 -> X3 -> X4")
```

## Response of single variable to two periodic forcings

```@docs
CausalityTools.Systems.verdes(;u₀ = rand(3),
    ωy = 315, ωz = 80,
    σx = 0.0, σy = 0.0, σz = 0.0)
```

```@example s
# Create an instance of the discrete dynamical system and
# produce an orbit of 100 points.
ds = CausalityTools.Systems.logistic4()
orbit = trajectory(ds, 100)

# Plot the orbit
t = 1:size(orbit, 1)
plot(t, orbit[:, 1], label = "X")
plot!(t, orbit[:, 2], label = "Y")
plot!(t, orbit[:, 3], label = "Z")
xlabel!("Time step")
ylabel!("Value")
title!("Periodic forcings from Y->X and Z->X")
```
