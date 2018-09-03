# Fourier phase surrogates


```@setup s
using CausalityTools
using DynamicalSystems
using Plots
```

Generating random phase surrogates is done as follows.

```@example s
# Generate a dynamical system, create an orbit and extract a time series.
s = CausalityTools.Systems.logistic3()
orbit = trajectory(s, 150)

# Extract our time series. We'll use the second component of the orbit.
x₂ = orbit[:, 2]

# Construct a random shuffle surrogate.
s₂ = randomphases(x₂)

# Plot the time series along with its random shuffle surrogate
timesteps = 1:size(orbit, 1)
plot(timesteps, x₂, label = "x2")
plot!(timesteps, s₂, label = "Random phase surrogate realization of x2")
```
