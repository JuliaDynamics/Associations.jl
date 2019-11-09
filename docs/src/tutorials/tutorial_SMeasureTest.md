# [Tutorial: SMeasureTest](@id tutorial_smeasuretest)

Create an orbit of the built-in unidirectionally coupled `henon2` map system, and 
a pair of random time series.

```julia
npts, Ttr = 5000, 500
x, y = columns(trajectory(henon2(c_xy = 1.0), npts, Ttr = Ttr));
xr, yr = rand(npts), rand(npts)
```

Initialise test, specifying embedding dimension, emebdding lag and number of 
nearest neighbors.

```julia
ks = 2:10
test = SMeasureTest(m = 4, τ = 1, K = ks)
```

Compute S-measure statistic separately in both directions, both for the 
random time series, and for the Henon maps. The test will return a vector 
of length `ks`.

```julia
Ss_r_xy = causality(xr, yr, test)
Ss_r_yx = causality(yr, xr, test)
Ss_henon_xy = causality(x, y, test)
Ss_henon_yx = causality(y, x, test);
```

Plot the results.

```julia
plot(xlabel = "# nearest neighbors (k)", ylabel = "S", ylims = (-0.05, 1.05))
plot!(Ks, Ss_r_xy,  label = "random uncoupled system (x -> y)", marker = stroke(2), c = :black)
plot!(Ks, Ss_r_yx,  label = "random uncoupled system (y -> x)", marker = stroke(2), c = :red)
plot!(Ks, Ss_henon_xy, marker = stroke(2), label = "henon unidir (x -> y)")
plot!(Ks, Ss_henon_yx, marker = stroke(2), label = "henon unidir (y -> x)")
```

![](figs/Example_Smeasure_random_and_henon.svg)