# Pre-defined coupled dynamical systems

CausalityTools provides a range of coupled dynamical systems, previously used in the literature to
study the behaviour of various causality statistics.

## Discrete systems 

### Lattice of unidirectionally coupled maps 

```@docs 
latticeunidir
```

```@example 
using CausalityTools, Plots, LaTeXStrings
sys = latticeunidir(50)
npts = 100
tr = trajectory(sys, npts)

# Plot time series for the first two variables
x1, x2 = tr[:, 1], tr[:, 2]
plot(xlabel = "Time", ylabel = "Value")
plot!(x1, label = L"x_1", c = :black)
plot!(x2, label = L"x_2", c = :red)
```

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

# Plot time series for the first component of each subsystem
x, y = orbit[:, 1], orbit[:, 4]
p_ts = plot(xlabel = "Time step", ylabel = "Value", bg_legend = :transparent)
plot!(x .- mean(x), label = L"x_1", c = :black)
plot!(y .- mean(y), label = L"y_1", c = :red)
```