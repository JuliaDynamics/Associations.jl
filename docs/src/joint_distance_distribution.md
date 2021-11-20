
# [Joint distance distribution](@id joint_distance_distribution_overview)

```@docs
jdd(::Any, ::Any)
```

For the joint distance distribution to indicate a causal influence, it must be significantly 
skewed towards positive values.

Providing the `OneSampleTTest` type as the first 
argument to `jdd` yields a one sample t-test on the joint distance distribution. From this test, you can extract p-values and obtain 
confidence intervals like in [HypothesisTests.jl](https://github.com/JuliaStats/HypothesisTests.jl) as usual.

```@docs
jdd(::Type{OneSampleTTest}, ::Any, ::Any)
```

## Example: two coupled Lorenz attractors

The joint distance distribution was introduced in Amigo et al. (2018)[^Amigo2018]. Here, we'll attempt to reproduce 
figures 1a and 1b from their paper, but first we'll do a couple of simple examples.

Let's start by generating some time series. We'll use the `lorenz_lorenz_bidir` example system that ships with CausalityTools, which implements their bidirectionally coupled Lorenz systems. 

```@example jdd
using CausalityTools, Plots, DynamicalSystems

# Time series of length 10000, sampled every 0.1 time steps.
ΔT, npts, Ttr = 0.1, 1000, 100
T = npts * ΔT

# Let there be unidirectional coupling from x to y
cxy, cyx = 1.5, 0.0
sys = lorenz_lorenz_bidir(c_xy = cxy, c_yx = cyx)
orbit = trajectory(sys, T, ΔT = ΔT, Ttr = 500)
x1, x2, x3, y1, y2, y3 = columns(orbit)
plot(xlabel = "Time step", ylabel = "Value", size = (800, 220))
plot!(x1, label = "x1")
plot!(y1, label = "y1")
savefig("lorenz_lorenz_timeseries.svg") # hide
```

![](lorenz_lorenz_timeseries.svg)

Let's compute the joint distance distribution from `x` to `y`. For this example, we'll use an embedding dimension 
of 10 and embedding lag of 10. In real applications, these parameters should be optimized.

```@example jdd
jxy = jdd(x1, y1, D = 10, τ = 10)
```

If `jxy` is biased towards positive values, then in the framework of the joint distance distribution, the variables are coupled and there exists a continuous functional dependence ``x = \phi(y)``. In the context of causal inference, this is is the same as saying that `x` drives `y`. Likewise, if `jxy` is *not* biased towards
positive values, then there exists no such relationship, and the variables are not coupled.

To formally test whether the distribution `jxy` is skewed torwards positive values, we can use a one-sample t-test to test the null hypothesis `mean(jxy) == 0`, and compute the right-sided p-value for the test. If `p < 0.05` for the test, we take that as rejection of the null hypothesis, and accept the alternative hypothesis that `mean(jxy)` is skewed towards positive values.

```@example jdd
pvalue(jdd(OneSampleTTest, x1, y1, D = 10, τ = 10), tail = :right)
```

Similarly, we can also see if the joint distance distribution test provides evidence of a coupling from `y` to `x`:

```@example jdd
pvalue(jdd(OneSampleTTest, y1, x1, D = 10, τ = 10), tail = :right)
```

If the test works perfectly, then we expect `p < 0.05` when computing the joint distance distribution from `x` to `y` and `p >= 0.05` when computing the distribution from `y` to `x`. However, for a single finite time series realization, this might not necessarily always be the case.  

## Reproducing Amigo et al. (2018)

Here, we attempt reproduce figures 1a and 1b from Amigo et al., which computes p-values for a range of coupling strengths in both directions. We'll use a coarser coupling resolution and shorter time series to limit computation time.

```@example jdd2
using CausalityTools, Plots, DynamicalSystems

cxys = 0.0:0.25:2.0
cyxs = 0.0:0.25:2.0
# We'll use time series 2000 points long, to limit computation time
ΔT, npts, Ttr = 0.1, 200, 100
T = npts * ΔT

pvals_xy = zeros(length(cyxs), length(cxys))
pvals_yx = zeros(length(cyxs), length(cxys))

for (i, cxy) in enumerate(cxys)
    for (j, cyx) in enumerate(cyxs)
        sys = lorenz_lorenz_bidir(c_xy = cxy, c_yx = cyx)
        orbit = trajectory(sys, T, ΔT = ΔT, Ttr = Ttr)
        x1, x2, x3, y1, y2, y3 = columns(orbit)

        # The original paper uses an embedding dimension of 10 and embedding lag of 10
        pvals_xy[j, i] = pvalue(jdd(OneSampleTTest, x1, y1, D = 10, τ = 10), tail = :right)
        pvals_yx[j, i] = pvalue(jdd(OneSampleTTest, y1, x1, D = 10, τ = 10), tail = :right)
    end
end
```

Plotting the p-values for the tests in the direction `x` to `y`, and in the direction `y` to `x`, as a function
of coupling strengths in both directions, we get the following:

```@example jdd2
z = [0, 0.0001, 0.001, 0.01, 0.05, 1.0]
c = cgrad(:magma, categorical = true)

p_xy = plot(xlabel = "Coupling from x to y", ylabel = "Coupling from y to x")
yticks!((1:length(cyxs), string.(cyxs)))
xticks!((1:length(cxys), string.(cxys)))
heatmap!(pvals_xy, c = c, logscale = true)

p_yx = plot(xlabel = "Coupling from x to y", ylabel = "")
yticks!((1:length(cyxs), string.(cyxs)))
xticks!((1:length(cxys), string.(cxys)))
heatmap!(pvals_yx, c = c, logscale = true)

plot(p_xy, p_yx, layout = grid(1, 2), size = (900, 300))
savefig("jdd_heatmap.svg") # hide
```

![](jdd_heatmap.svg)

The colorbar indicates the p-value. For sufficiently high coupling strengths, the test successfully identifies the couplings both from `x` to `y` and vice versa. For lower coupling strengths, we cannot reject the null hypothesis. 

Note that using longer time series and more carefully tuned embeddings would likely lead to better performance.


[^Amigo2018]: Amigó, José M., and Yoshito Hirata. "Detecting directional couplings from multivariate flows by the joint distance distribution." Chaos: An Interdisciplinary Journal of Nonlinear Science 28.7 (2018): 075302.
