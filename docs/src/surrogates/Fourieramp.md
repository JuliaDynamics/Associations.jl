# Fourier amplitude surrogates

```@setup s
using CausalityTools
using DynamicalSystems
using Plots
```

Generating random amplitude Fourier surrogates is done as follows.

```@example s
# Generate a dynamical system, create an orbit and extract a time series.
s = CausalityTools.Systems.logistic3()
orbit = trajectory(s, 150)

# Extract our time series. We'll use the first component of the orbit.
x₁ = orbit[:, 1]

# Construct a random shuffle surrogate.
s₁ = randomamplitudes(x₁)

# Plot the time series along with its random shuffle surrogate
timesteps = 1:size(orbit, 1)
plot(timesteps, x₁, label = "x1")
plot!(timesteps, s₁, label = "Random amplitude surrogate realization of x1")
xlabel!("Timestep"); ylabel!("Value")
savefig("surr_randomamplitude-plot.svg"); nothing #hide
```

![](surr_randomamplitude-plot.svg)
