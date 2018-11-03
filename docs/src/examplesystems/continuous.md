
# Continuous systems

Where stated, we give instructions on what solver and which time steps the
original authors used. You can pass these along as keyword arguments to the
`trajectory` function from `DynamicalSystems.jl`, which we use to generate
orbits.

```@setup s
using CausalityTools
using Plots
using DynamicalSystems
```

## Triple Lorenz system

```@docs
CausalityTools.Systems.lorenz_triple(;uᵢ=rand(9),
        σ₁ = 10.0, σ₂ = 10.0, σ₃ = 10.0,
        ρ₁ = 28.0, ρ₂ = 28.0, ρ₃ = 28.0,
        β₁ = 8/3, β₂ = 8/3, β₃ = 8.3,
        c₁ = 1.0, c₂ = 1.0)
```

```@example s
# Create an instance of the continuous dynamical system and integrate
# it up to time `t = 500` in steps of `dt = 0.01`. This given an
# orbit of 50000 points. We sample every fifth point of the orbit.
ds = CausalityTools.Systems.lorenz_triple()
orbit = trajectory(ds, 200, dt = 0.01)[1:5:end, :]

# Plot the orbit
t = 1:size(orbit, 1)
plot(t, orbit[:, 2], label = "y-component of X1")
plot!(t, orbit[:, 5], label = "y-component of X2")
plot!(t, orbit[:, 8], label = "y-component of X3")
xlabel!("Time step"); ylabel!("Value")
title!("Three coupled Lorenz systems, X1 -> X2 and X2 -> X3")
savefig("triplelorenz-plot.svg"); nothing #hide
```

![](triplelorenz-plot.svg)


## Rössler-Lorenz system

```@docs
rossler_lorenz(;u₀ = rand(6), a₁ = -6, a₂ = 6,
							a₃ = 2.0, b₁ = 10, b₂ = 28, b₃ = 8/3, c = 2)
```


```@example s
# Create an instance of the continuous dynamical system and integrate
# it up to time `t = 300` in steps of `dt = 0.01`. We sample every 6th
# point of the orbit.
ds = CausalityTools.Systems.rossler_lorenz()
orbit = trajectory(ds, 300, dt = 0.01)[1:6:end, :]

# Plot the orbit
t = 1:size(orbit, 1)
plot(t, orbit[:, 3], label = "Z1")
plot!(t, orbit[:, 6], label = "Z2")
xlabel!("Time step"); ylabel!("Value")
title!("Rössler-Lorenz system, Z1 -> Z2")
savefig("rosslerlorenz-plot.svg"); nothing #hide
```

![](rosslerlorenz-plot.svg)



## Rössler-Rössler system


```@docs
  rossler_rossler(;uᵢ = rand(6), a = 0.1, b = 0.1, c = 14.0, ϵ₁ = 0.0, ϵ₂ = 0.0, ω₁ = 1 + 0.015, ω₂ = 1 - 0.015)
```

Let's create a system of coupled Rössler attractors with different frequencies.

```@example s
# Create an instance of the continuous dynamical system and integrate
# it up to time `t = 200` in steps of `dt = 0.01`. We sample every 3rd
# point of the orbit, and discard points corresponding to times up to t = 10.
ds = CausalityTools.Systems.rossler_rossler(ω₁ = 1.5, ω₂ = 0.5)
orbit = trajectory(ds, 200, dt = 0.01)[10001:3:end, :]

# Plot the orbit
t = 1:size(orbit, 1)
plot(t, orbit[:, 1], label = "X1")
plot!(t, orbit[:, 4], label = "X2")
xlabel!("Time step"); ylabel!("Value")
title!("Rössler-Rössler system, X1 -> X2")
savefig("rosslerrossler-plot.svg"); nothing #hide
```

![](rosslerrossler-plot.svg)
