# Pre-defined coupled dynamical systems

CausalityTools provides a range of coupled dynamical systems, previously used in the literature to
study the behaviour of various causality statistics.

## Continuous systems 

### Lorenz, diffusively coupled

```@docs
lorenzdiffusive
```

```@example
using CausalityTools, DynamicalSystemsBase, Statistics, Plots, LaTeXStrings; pyplot()
ld = lorenzdiffusive()
npts, nt, Δt = 3000, 10000, 0.03
orbit = trajectory(ld, Δt*(npts-1), dt = Δt, Ttr = nt*Δt)
x, y = orbit[:, 1], orbit[:, 4]

p_ts = plot(xlabel = "Time step", ylabel = "Value", bg_legend = :transparent)
plot!(x .- mean(x), label = L"x_1", c = :black)
plot!(y .- mean(y), label = L"y_1", c = :red)
```