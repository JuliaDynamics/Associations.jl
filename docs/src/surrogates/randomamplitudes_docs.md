# Surrogate methods

## Random amplitude Fourier surrogates

```@docs
randomamplitudes
```

### Example

```@setup randomamplitudes
using CausalityTools, Plots
```

```@example randomamplitudes
# Generate a dynamical system, create an orbit and extract a time series from
# the first component.
s = CausalityTools.Systems.logistic3()
orbit = trajectory(s, 150)
x = orbit[:, 1]

# Compare original time series and random amplitude surrogate
p1 = plot(x, label = "x", lc = :black)
p2 = plot(randomamplitudes(x), label = "randomamplitudes(x)",
            xlabel = "Time step")
plot(p1, p2, layout = (2, 1))
ylabel!("Value")
savefig("surr_randomamplitudes.svg"); nothing #hide
```

![](surr_randomamplitudes.svg)
