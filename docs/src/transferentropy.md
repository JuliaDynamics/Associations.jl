# Transfer entropy

The transfer entropy API is made up of the following functions and types,
which are listed below:

- [`transferentropy`](@ref).
- [`TransferEntropy`](@ref), and its subtypes.
- [`EmbeddingTE`](@ref), which exists to provide embedding instructions to
    subtypes of [`TransferEntropy`](@ref).
- [`TransferEntropyEstimator`](@ref), and its subtypes.

## API

```@docs
transferentropy
EmbeddingTE
```

## Definitions

```@docs
TransferEntropy
TEShannon
TERenyiJizba
```

## Estimators

```@docs
TransferEntropyEstimator
Zhu1
Lindner
```

## Examples

Let's try to reproduce the results from Schreiber's original paper[^Schreiber2000] where
he introduced the transfer entropy. We'll use the [`ValueHistogram`](@ref) estimator,
which is visitation frequency based and computes entropies by counting visits of the
system's orbit to discrete portions of its reconstructed state space.

```@example example_te_schreiber
using CausalityTools
using DynamicalSystemsBase
using CairoMakie
using Statistics
using Random; Random.seed!(12234);

function ulam_system(dx, x, p, t)
    f(x) = 2 - x^2
    ε = p[1]
    dx[1] = f(ε*x[length(dx)] + (1-ε)*x[1])
    for i in 2:length(dx)
        dx[i] = f(ε*x[i-1] + (1-ε)*x[i])
    end
end

ds = DiscreteDynamicalSystem(ulam_system, rand(100) .- 0.5, [0.04])
trajectory(ds, 1000; Ttr = 1000);

εs = 0.02:0.02:1.0
base = 2
te_x1x2 = zeros(length(εs)); te_x2x1 = zeros(length(εs))
# Guess an appropriate bin width of 0.2 for the histogram
est = ValueHistogram(0.2)

for (i, ε) in enumerate(εs)
    set_parameter!(ds, 1, ε)
    tr = trajectory(ds, 2000; Ttr = 5000)
    X1 = tr[:, 1]; X2 = tr[:, 2]
    @assert !any(isnan, X1)
    @assert !any(isnan, X2)
    te_x1x2[i] = transferentropy(TEShannon(; base), est, X1, X2)
    te_x2x1[i] = transferentropy(TEShannon(; base), est, X2, X1)
end

fig = with_theme(theme_minimal(), markersize = 2) do
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "epsilon", ylabel = "Transfer entropy (bits)")
    scatterlines!(ax, εs, te_x1x2, label = "X1 to X2", color = :black, lw = 1.5)
    scatterlines!(ax, εs, te_x2x1, label = "X2 to X1", color = :red, lw = 1.5)
    axislegend(ax, position = :lt)
    return fig
end
fig
```

As expected, transfer entropy from `X1` to `X2` is higher than from `X2` to `X1` across parameter values for `ε`. But, by our definition of the ulam system, dynamical coupling only occurs from `X1` to `X2`. The results, however, show nonzero transfer entropy in both directions. What does this mean?

Computing transfer entropy from finite time series introduces bias, and so does any particular choice of entropy estimator used to calculate it. To determine whether a transfer entropy estimate should be trusted, we can employ surrogate testing. We'll generate surrogate using
[TimeseriesSurrogates.jl](https://github.com/JuliaDynamics/TimeseriesSurrogates.jl).

In the example below, we continue with the same time series generated above. However, at each value of `ε`, we also compute transfer entropy for `nsurr = 50` different randomly shuffled (permuted) versions of the source process. If the original transfer entropy exceeds that of some percentile the transfer entropy estimates of the surrogate ensemble, we will take that as "significant" transfer entropy.

```@example example_te_schreiber
nsurr = 25 # in real applications, you should use more surrogates
base = 2
te_x1x2 = zeros(length(εs)); te_x2x1 = zeros(length(εs))
te_x1x2_surr = zeros(length(εs), nsurr); te_x2x1_surr = zeros(length(εs), nsurr)
est = ValueHistogram(0.2) # use same bin-width as before

for (i, ε) in enumerate(εs)
    set_parameter!(ds, 1, ε)
    tr = trajectory(ds, 500; Ttr = 5000)
    X1 = tr[:, 1]; X2 = tr[:, 2]
    @assert !any(isnan, X1)
    @assert !any(isnan, X2)
    te_x1x2[i] = transferentropy(TEShannon(; base), est, X1, X2)
    te_x2x1[i] = transferentropy(TEShannon(; base), est, X2, X1)
    s1 = surrogenerator(X1, RandomShuffle()); s2 = surrogenerator(X2, RandomShuffle())

    for j = 1:nsurr
        te_x1x2_surr[i, j] =  transferentropy(TEShannon(; base), est, s1(), X2)
        te_x2x1_surr[i, j] =  transferentropy(TEShannon(; base), est, s2(), X1)
    end
end

# Compute 95th percentiles of the surrogates for each ε
qs_x1x2 = [quantile(te_x1x2_surr[i, :], 0.95) for i = 1:length(εs)]
qs_x2x1 = [quantile(te_x2x1_surr[i, :], 0.95) for i = 1:length(εs)]

fig = with_theme(theme_minimal(), markersize = 2) do
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "epsilon", ylabel = "Transfer entropy (bits)")
    scatterlines!(ax, εs, te_x1x2, label = "X1 to X2", color = :black, lw = 1.5)
    scatterlines!(ax, εs, qs_x1x2, color = :black, linestyle = :dot, lw = 1.5)
    scatterlines!(ax, εs, te_x2x1, label = "X2 to X1", color = :red)
    scatterlines!(ax, εs, qs_x2x1, color = :red, linestyle = :dot)
    axislegend(ax, position = :lt)
    return fig
end
fig
```

The plot above shows the original transfer entropies (solid lines) and the 95th percentile transfer entropies of the surrogate ensembles (dotted lines). As expected, using the surrogate test, the transfer entropies from `X1` to `X2` are mostly significant (solid black line is above dashed black line). The transfer entropies from `X2` to `X1`, on the other hand, are mostly not significant (red solid line is below red dotted line).

[^Schreiber2000]:
    Schreiber, Thomas. "Measuring information transfer." Physical review letters 85.2
    (2000): 461.

### Estimator comparison

Let's reproduce Figure 4 from Zhu et al (2015)[^Zhu2015], where they test some
dedicated transfer entropy estimators on a bivariate autoregressive system.
We will test

- The [`Lindner`](@ref) and [`Zhu1`](@ref) dedicated transfer entropy estimators,
    which try to eliminate bias.
- The [`KSG1`](@ref) estimator, which computes TE naively as a sum of mutual information
    terms (without guaranteed cancellation of biases for the total sum).
- The [`Kraskov`](@ref) estimator, which computes TE naively as a sum of entropy 
    terms (without guaranteed cancellation of biases for the total sum).

[^Zhu2015]:
    Zhu, J., Bellanger, J. J., Shu, H., & Le Bouquin Jeannès, R. (2015). Contribution to transfer entropy estimation via the k-nearest-neighbors approach. Entropy, 17(6), 4173-4201.

```@example
using CausalityTools
using CairoMakie
using Statistics
using Distributions: Normal

function model2(n::Int)
    𝒩x = Normal(0, 0.1)
    𝒩y = Normal(0, 0.1)
    x = zeros(n+2)
    y = zeros(n+2)
    x[1] = rand(𝒩x)
    x[2] = rand(𝒩x)
    y[1] = rand(𝒩y)
    y[2] = rand(𝒩y)

    for i = 3:n+2
        x[i] = 0.45*sqrt(2)*x[i-1] - 0.9*x[i-2] - 0.6*y[i-2] + rand(𝒩x)
        y[i] = 0.6*x[i-2] - 0.175*sqrt(2)*y[i-1] + 0.55*sqrt(2)*y[i-2] + rand(𝒩y)
    end
    return x[3:end], y[3:end]
end
te_true = 0.42 # eyeball the theoretical value from their Figure 4.

m = TEShannon(embedding = EmbeddingTE(dT = 2, dS = 2), base = ℯ)
estimators = [Zhu1(k = 8), Lindner(k = 8), KSG1(k = 8), Kraskov(k = 8)]
Ls = [floor(Int, 2^i) for i in 8.0:0.5:11]
nreps = 8
tes_xy = [[zeros(nreps) for i = 1:length(Ls)] for e in estimators]
tes_yx = [[zeros(nreps) for i = 1:length(Ls)] for e in estimators]
for (k, est) in enumerate(estimators)
    for (i, L) in enumerate(Ls)
        for j = 1:nreps
            x, y = model2(L);
            tes_xy[k][i][j] = transferentropy(m, est, x, y)
            tes_yx[k][i][j] = transferentropy(m, est, y, x)
        end
    end
end

ymin = minimum(map(x -> minimum(Iterators.flatten(Iterators.flatten(x))), (tes_xy, tes_yx)))
estimator_names = ["Zhu1", "Lindner", "KSG1", "Kraskov"]
ls = [:dash, :dot, :dash, :dot]
mr = [:rect, :hexagon, :xcross, :pentagon]

fig = Figure(resolution = (800, 350))
ax_xy = Axis(fig[1,1], xlabel = "Signal length", ylabel = "TE (nats)", title = "x → y")
ax_yx = Axis(fig[1,2], xlabel = "Signal length", ylabel = "TE (nats)", title = "y → x")
for (k, e) in enumerate(estimators)
    label = estimator_names[k]
    marker = mr[k]
    scatterlines!(ax_xy, Ls, mean.(tes_xy[k]); label, marker)
    scatterlines!(ax_yx, Ls, mean.(tes_yx[k]); label, marker)
    hlines!(ax_xy, [te_true]; xmin = 0.0, xmax = 1.0, linestyle = :dash, color = :black) 
    hlines!(ax_yx, [te_true]; xmin = 0.0, xmax = 1.0, linestyle = :dash, color = :black)
    linkaxes!(ax_xy, ax_yx)
end
axislegend(ax_xy, position = :rb)

fig
```

## Convenience

### Symbolic transfer entropy

```@docs
SymbolicTransferEntropy
```

### Phase/amplitude transfer entropy

```@docs
Hilbert
Phase
Amplitude
```
