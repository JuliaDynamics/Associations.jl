# [S-measure](@id Smeasure_overview)

```@docs
s_measure
```
## Example: random data vs Henon maps

Here, we'll compute the S-measure between random time series (uniform noise), between
time series of a dynamical system (coupled Henon maps).

We start by generating the time series.

```@example smeasure_random_henon
using CausalityTools, DynamicalSystems, Plots, Distributions
npts, Ttr = 10000, 5000
x, y = columns(trajectory(henon2(c_xy = 0.87), npts - 1, Ttr = Ttr))
xr, yr = rand(Uniform(-1, 1), npts), rand(Uniform(-1, 1), npts)
[x y xr yr]
```

Let's plot the time series.

```@example smeasure_random_henon
p_det = Plots.plot(xlabel = "", ylabel = "Value", title = "Coupled Henon maps")
Plots.plot!(x[1:100], label = "x", marker = stroke(:black), c = :black)
Plots.plot!(y[1:100], label = "y", marker = stroke(:red), c = :red)
p_rand = Plots.plot(xlabel = "Time step", ylabel = "Value", title = "Random time series")
Plots.plot!(xr[1:100], label = "xr", c = :blue)
Plots.plot!(yr[1:100], label = "yr", c = :purple)

Plots.plot(p_det, p_rand, layout = Plots.grid(2, 1), size = (382*2, 400), legend = :bottomright, 
    tickfont = font(13), guidefont = font(13), legendfont = font(13))
savefig("henon_random_timeseries.svg"); nothing # hide
```

![](henon_random_timeseries.svg)

Now we compute the S-measure between the random time series, both from `x` to `y` and from `y` to `x`.
We'll also do the same for the Henon maps. 

The test parameters are embedding dimensions (`dx` for the source and `dy` for the target), the embedding lags (`τx` for the source and `τy` for the target), and the number of nearest neighbors `K`. We'll compute the test with fixed embedding parameters, but a varying number of nearest neighbors (`ks = 2:10`).

For the sake of demonstration, we'll use 4-dimensional embedding with embedding lag 3 for the source, and a 5-dimensional embedding with embedding lag 1 for the target. For a real use case, these embedding parameters should be 
chosen more carefully.

```@example smeasure_random_henon
ks = 2:8
# Compute the s-measures for different values of k
ss_r_xy = [s_measure(x, y, dx = 4, τx = 3, dy = 5, τy = 1, K = k) for k in ks]
ss_r_yx = [s_measure(yr, xr, dx = 4, τx = 3, dy = 5, τy = 1, K = k) for k in ks]
ss_henon_xy = [s_measure(x, y, dx = 4, τx = 3, dy = 5, τy = 1, K = k) for k in ks]
ss_henon_yx = [s_measure(y, x, dx = 4, τx = 3, dy = 5, τy = 1, K = k) for k in ks]
[ss_r_xy ss_r_yx ss_henon_xy ss_henon_yx]

Plots.plot(xlabel = "# nearest neighbors (k)", ylabel = "S", ylims = (-0.05, 1.05))
Plots.plot!(ks, ss_r_xy,  label = "random uncoupled system (x -> y)", marker = stroke(2), c = :black)
Plots.plot!(ks, ss_r_yx,  label = "random uncoupled system (y -> x)", marker = stroke(2), c = :red)
Plots.plot!(ks, ss_henon_xy, marker = stroke(2), label = "henon unidir (x -> y)")
Plots.plot!(ks, ss_henon_yx, marker = stroke(2), label = "henon unidir (y -> x)")
savefig("smeasure_random_henon.svg"); nothing # hide
```

![](smeasure_random_henon.svg)

For uncoupled time series, we expect the value of $S$ to be close to zero. For strongly coupled time series, the value of $S$ should be nonzero and approaching 1. This is exactly what we get: for time random time series, the value of $S$ is close to zero and for the Henon maps, it's clearly non-zero.

Note that the actual dynamical coupling in the Henon system is unidirectional from `x` to `y`. The results (positive $S$ in both directions), however, indicate that the coupling is bidirectional, with the coupling being stronger in one direction. This disagreement between results and ground truth highlights the importance of employing *causal hypothesis testing*. For this, we could use [TimeseriesSurrogates.jl](https://github.com/JuliaDynamics/TimeseriesSurrogates.jl) to 
generate surrogate time series. An in-depth example on how to use surrogate testing can be found in the [transfer entropy](@ref transferentropy) example.
