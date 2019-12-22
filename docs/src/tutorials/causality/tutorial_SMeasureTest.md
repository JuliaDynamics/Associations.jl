# [Tutorial: SMeasureTest](@id tutorial_smeasuretest)

Create an orbit of the built-in unidirectionally coupled `henon2` map system, and 
a pair of random time series.

```@example SMeasureTest_henon2_rand
using CausalityTools, DynamicalSystemsBase, Plots
npts, Ttr = 5000, 500
x, y = columns(trajectory(henon2(c_xy = 1.0), npts, Ttr = Ttr));
xr, yr = rand(npts), rand(npts); nothing # hide
```

Initialise test, specifying embedding dimension, emebdding lag and number of 
nearest neighbors.

```@example SMeasureTest_henon2_rand
ks = 2:10
test = SMeasureTest(m = 4, τ = 1, K = ks)
```

Compute S-measure statistic separately in both directions, both for the 
random time series, and for the Henon maps. The test will return a vector 
of length `ks`.

```@example SMeasureTest_henon2_rand
Ss_r_xy = causality(xr, yr, test)
Ss_r_yx = causality(yr, xr, test)
Ss_henon_xy = causality(x, y, test)
Ss_henon_yx = causality(y, x, test);
nothing # hide
```

Plot the results.

```@example SMeasureTest_henon2_rand
plot(xlabel = "# nearest neighbors (k)", ylabel = "S", ylims = (-0.05, 1.05))
plot!(ks, Ss_r_xy,  label = "random uncoupled system (x -> y)", marker = stroke(2), c = :black)
plot!(ks, Ss_r_yx,  label = "random uncoupled system (y -> x)", marker = stroke(2), c = :red)
plot!(ks, Ss_henon_xy, marker = stroke(2), label = "henon unidir (x -> y)")
plot!(ks, Ss_henon_yx, marker = stroke(2), label = "henon unidir (y -> x)")
savefig("figs/SMeasure_random_plus_henon.svg"); nothing # hide
```

![](figs/SMeasure_random_plus_henon.svg)


For uncoupled time series, we expect the value of $S$ to be close to zero. For strongly coupled time series, the value of $S$ should be nonzero and approaching 1. This is exactly what we get: for time random time series, the value of $S$ is close to zero and for the Henon maps, it's clearly non-zero.