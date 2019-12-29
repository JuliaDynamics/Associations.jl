# [`SMeasureTest`](@id SMeasureTest_overview)

```@docs
SMeasureTest
```

## Example

### Time series

First, create an orbit of the built-in unidirectionally coupled [`henon2`](@ref) map system,
and a set of random time series for comparison.

```@example SMeasureTest_henon2_rand
using CausalityTools, DynamicalSystems, Plots, Distributions

npts, Ttr = 5000, 500
x, y = columns(trajectory(henon2(c_xy = 1.0), npts - 1, Ttr = Ttr))
xr, yr = rand(Uniform(-1, 1), npts), rand(Uniform(-1, 1), npts)
[x y xr yr]
```

Let's plot the first 100 points of each time series.

```@example SMeasureTest_henon2_rand
p_det = plot(xlabel = "", ylabel = "Value", title = "Coupled Henon maps")
plot!(x[1:100], label = "x", marker = stroke(:black), c = :black)
plot!(y[1:100], label = "y", marker = stroke(:red), c = :red)
p_rand = plot(xlabel = "Time step", ylabel = "Value", title = "Random time series")
plot!(xr[1:100], label = "xr", c = :blue)
plot!(yr[1:100], label = "yr", c = :purple, ls = :dash)

plot(p_det, p_rand, layout = grid(2, 1), size = (382*2, 400), legend = :bottomright,
    tickfont = font(13), guidefont = font(13), legendfont = font(13))
```

### Test setup

Create a [`SMeasureTest`](@ref), which contains the test parameters. These are 
the embedding dimension `m`, the embedding lag `τ`, and the number of nearest neighbors `K`. 
We'll compute the test with fixed embedding, but a varying number of nearest neighbors 
`ks = 2:10`.

```@example SMeasureTest_henon2_rand
ks = 2:10
test = SMeasureTest(m = 4, τ = 1, K = ks)
```

### Analysis

Compute S-measure statistic separately in both directions, both for the 
random time series, and for the Henon maps. The test will return a vector 
of length `ks`.

```@example SMeasureTest_henon2_rand
Ss_r_xy = causality(xr, yr, test)
Ss_r_yx = causality(yr, xr, test)
Ss_henon_xy = causality(x, y, test)
Ss_henon_yx = causality(y, x, test);
[Ss_r_xy Ss_r_yx Ss_henon_xy Ss_henon_yx]
```

The first two columns contain the s-measures computed for our deterministic Henon map 
time series, and the two last columns contain the s-measures computed for the stochastic
time series.

### Results

Let's plot the results.

```@example SMeasureTest_henon2_rand
plot(xlabel = "# nearest neighbors (k)", ylabel = "S", ylims = (-0.05, 1.05))
plot!(ks, Ss_r_xy,  label = "random uncoupled system (x -> y)", marker = stroke(2), c = :black)
plot!(ks, Ss_r_yx,  label = "random uncoupled system (y -> x)", marker = stroke(2), c = :red)
plot!(ks, Ss_henon_xy, marker = stroke(2), label = "henon unidir (x -> y)")
plot!(ks, Ss_henon_yx, marker = stroke(2), label = "henon unidir (y -> x)")
```

### Discussion

For uncoupled time series, we expect the value of $S$ to be close to zero. For strongly coupled time series, the value of $S$ should be nonzero and approaching 1. This is exactly what we get: for time random time series, the value of $S$ is close to zero and for the Henon maps, it's clearly non-zero.
