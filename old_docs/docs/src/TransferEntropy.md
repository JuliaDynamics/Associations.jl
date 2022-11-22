# [Transfer entropy](@ref transferentropy)

## Traditional 

The following `transferentropy` function computes transfer entropy "manually", that is,
in addition to specifying an estimator, you have to specify embedding parameters.

```@docs
transferentropy
```

## Automated variable selection

The `bbnue` function optimizes embedding parameters using an iterative procedure for 
variable selection, and performs null hypothesis testing as part of that procedure.

```@docs
bbnue
```

## Example: Reproducing Schreiber (2000)

Let's try to reproduce the results from Schreiber's original paper[^Schreiber2000] on transfer entropy. We'll use a 
visitation frequency estimator, which computes entropies by counting visits of the system's orbit to discrete portions 
of its reconstructed state space.

```@example schreiber1
using DynamicalSystems, CausalityTools, Plots, Random, StatsBase

Random.seed!(12234)

function ulam_system(dx, x, p, t)
    f(x) = 2 - x^2
    ε = p[1]
    dx[1] = f(ε*x[length(dx)] + (1-ε)*x[1])
    for i in 2:length(dx)
        dx[i] = f(ε*x[i-1] + (1-ε)*x[i])
    end
end

ds = DiscreteDynamicalSystem(ulam_system, rand(100) .- 0.5, [0.04])
trajectory(ds, 1000; Ttr = 1000)

εs = 0.02:0.02:1.0
te_x1x2 = zeros(length(εs)); te_x2x1 = zeros(length(εs))

for (i, ε) in enumerate(εs)
    set_parameter!(ds, 1, ε)
    tr = trajectory(ds, 2000; Ttr = 5000)
    X1 = tr[:, 1]; X2 = tr[:, 2]
    @assert !any(isnan, X1)
    @assert !any(isnan, X2)
    binning = RectangularBinning(0.2) # guess an appropriate bin width of 0.2
    te_x1x2[i] = transferentropy(X1, X2, VisitationFrequency(binning), base = 2)
    te_x2x1[i] = transferentropy(X2, X1, VisitationFrequency(binning), base = 2)
end

plot()
plot(εs, te_x1x2, label = "X1 to X2", c = :black, lw = 1.5)
plot!(εs, te_x2x1, label = "X2 to X1", c = :red)
xlabel!("epsilon")
ylabel!("Transfer entropy (bits)")
savefig("ulam-te.svg"); nothing # hide
```

![](ulam-te.svg)

As expected, transfer entropy from `X1` to `X2` is higher than from `X2` to `X1` across parameter values for `ε`. But, by our definition of the `ulam` system, dynamical coupling only occurs from `X1` to `X2`. The results, however, 
show nonzero transfer entropy in both directions. What does this mean? 

Computing transfer entropy from finite time series introduces bias, and so does any particular choice of entropy estimator used to calculate it. To determine whether a transfer entropy estimate should be trusted, we can employ
*surrogate testing*. We'll generate surrogate using [TimeseriesSurrogates.jl](https://github.com/JuliaDynamics/TimeseriesSurrogates.jl).

In the example below, we continue with the same time series generated above. However, at each value of `ε`, we also compute transfer entropy for `nsurr = 50` different randomly shuffled (permuted) versions of the source process. 
If the original transfer entropy exceeds that of some percentile the transfer entropy estimates of the surrogate ensemble, we will take that as "significant" transfer entropy.

```@example schreiber1
nsurr = 50
te_x1x2 = zeros(length(εs)); te_x2x1 = zeros(length(εs))
te_x1x2_surr = zeros(length(εs), nsurr); te_x2x1_surr = zeros(length(εs), nsurr)

for (i, ε) in enumerate(εs)
    set_parameter!(ds, 1, ε)
    tr = trajectory(ds, 1000; Ttr = 5000)
    X1 = tr[:, 1]; X2 = tr[:, 2]
    @assert !any(isnan, X1)
    @assert !any(isnan, X2)
    binning = RectangularBinning(0.2) # guess an appropriate bin width of 0.2
    est = VisitationFrequency(binning)
    te_x1x2[i] = transferentropy(X1, X2, est, base = 2)
    te_x2x1[i] = transferentropy(X2, X1, est, base = 2)
    s1 = surrogenerator(X1, RandomShuffle()); s2 = surrogenerator(X2, RandomShuffle())

    for j = 1:nsurr
        te_x1x2_surr[i, j] =  transferentropy(s1(), X2, est, base = 2)
        te_x2x1_surr[i, j] =  transferentropy(s2(), X1, est, base = 2)
    end
end

# Compute 95th percentiles of the surrogates for each ε
qs_x1x2 = [quantile(te_x1x2_surr[i, :], 0.95) for i = 1:length(εs)]
qs_x2x1 = [quantile(te_x2x1_surr[i, :], 0.95) for i = 1:length(εs)]

plot(xlabel = "epsilon", ylabel = "Transfer entropy (bits)", legend = :topleft)
plot!(εs, te_x1x2, label = "X1 to X2", c = :black, lw = 1.5)
plot!(εs, qs_x1x2, label = "", c = :black, ls = :dot, lw = 1.5)
plot!(εs, te_x2x1, label = "X2 to X1", c = :red)
plot!(εs, qs_x2x1, label = "", c = :red, ls = :dot)

savefig("ulam-te-surr.svg"); nothing # hide
```

![](ulam-te-surr.svg)


The plot above shows the original transfer entropies (solid lines) and the 95th percentile transfer entropies of the surrogate ensembles (dotted lines). As expected, using the surrogate test, the transfer entropies from `X1` to `X2` are mostly significant (solid black line is *above* dashed black line). The transfer entropies from `X2` to `X1`, on the other hand, are mostly not significant (red solid line is *below* red dotted line).

[^Schreiber2000](Schreiber, Thomas. "Measuring information transfer." Physical review letters 85.2 (2000): 461.)
