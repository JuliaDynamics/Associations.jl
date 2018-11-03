# Random shuffle surrogates

```@setup s
using CausalityTools
using DynamicalSystems
using Plots
```

Generating random shuffle surrogates is done as follows.

```@example s
# Generate a dynamical system, create an orbit and extract a time series.
s = CausalityTools.Systems.logistic4()
orbit = trajectory(s, 150)

# Extract our time series. We'll use the third component of the orbit.
x₃ = orbit[:, 3]

# Construct a random shuffle surrogate.
s₃ = randomshuffle(x₃)

# Plot the time series along with its random shuffle surrogate
timesteps = 1:size(orbit, 1)
plot(timesteps, x₃, label = "x3")
plot!(timesteps, s₃, label = "Random shuffle surrogate realization of x3")
xlabel!("Timestep"); ylabel!("Value")
savefig("surr_randomshuffle-plot.svg"); nothing #hide
```

![](surr_randomshuffle-plot.svg)
