# Discrete maps

```@setup s
using CausalityTools
using Plots
using DynamicalSystems
```

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
xlabel!("Time step"); ylabel!("Value")
title!("Two coupled logistic maps, forcing X -> Y")
savefig("logistic2-plot.svg"); nothing #hide
```

![](logistic2-plot.svg)

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
xlabel!("Time step"); ylabel!("Value")
title!("Three coupled logistic maps, Z -> X and Z -> Y")
savefig("logistic3-plot.svg"); nothing #hide
```

![](logistic3-plot.svg)





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
xlabel!("Time step"); ylabel!("Value")
title!("Periodic forcings from Y->X and Z->X")
savefig("verdes-plot.svg"); nothing #hide
```

![](verdes-plot.svg)
