# Surrogate methods

## Random phases Fourier surrogates

```@docs
randomphases
```

### Example

```@setup s
using CausalityTools
using DynamicalSystems
using Plots
```

Generating random phase surrogates is done as follows.

```@example s
# Generate a dynamical system, create an orbit and extract a time series from
# the second component.
s = CausalityTools.Systems.logistic3()
orbit = trajectory(s, 150)
x = orbit[:, 2]


# Plot the time series along with its random phase surrogate
timesteps = 1:size(orbit, 1)
p1 = plot(timesteps, x, label = "x")
p2 = plot(timesteps, randomphases(x), label = "randomphases(x)",
    xlabel = "Time step")
plot(p1, p2, layout = (2, 1))
ylabel!("Value")
savefig("surr_randomphases.svg"); nothing #hide
```

![](surr_randomphases.svg)

