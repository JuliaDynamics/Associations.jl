# [Generalized entropies](@id entropies)

## API

The following functions/types make up the entropy API:

- [`Entropy`](@ref), and all its subtypes.
- [`entropy`](@ref) / [`entropy!`](@ref)
- [`entropy_maximum`](@ref),
- [`entropy_normalized`](@ref),
- [`EntropyEstimator`](@ref), and all its subtypes.

```@docs
entropy
entropy!
entropy_maximum
entropy_normalized
```

## Discrete

Discrete entropies are defined by formulas defined over a set of [`Probabilities`](@ref), has a minimum value of 0, often have a defined maximum (see [`entropy_maximum`](@ref)), and can therefore be normalized (see [`entropy_normalized`](@ref)).

```@docs
Entropy
Renyi
Tsallis
Shannon
Curado
StretchedExponential
```

## Continuous/differential entropy

Continuous (differential) entropies are defined by an *integral*, and are 
related to but don't share all the same properties as discrete entropies.
For example, Shannon differential entropy may even be *negative* for 
some distributions.

Continuous entropies must be estimated using some form of "plug-in" estimator. For example, the Shannon differential entropy for a random variable $X$ with support $\mathcal{X}$ is defined as

$$h(x) = \mathbb{E}[-\log{(f(X))}] = -\int_{\mathcal{X}}f(x) \log f(x) dx.$$

There are several ways of estimating this integral from observed data,
using what is called "plug-in" estimators. A common plug-in estimator is the resubstitution estimator

$$\hat{H}(x) = -\frac{1}{N}\sum_{i=1}^N \log{(\hat{p}(X_i))},$$

where $\hat{p}$ is estimated using the samples $X_1, X_2, \ldots, X_N$, is a plug-in estimator for Shannon differential entropy.

Subtypes of [`EntropyEstimator`](@ref)s use various forms of plug-in estimators to estimate differential entropy. For example, [`Kraskov`](@ref) estimates *Shannon* differential entropy. [`LeonenkoProzantoSavani`](@ref), on the other hand, estimates both [`Shannon`](@ref), [`Renyi`](@ref) and 
[`Tsallis`](@ref) differential entropy.

!!! note "Plug-in estimators for differential entropy"
    When using [`entropy`](@ref) with a [`ProbabilityEstimator`](@ref), it is always the discrete entropy that is computed. When using [`entropy`](@ref) with an [`EntropyEstimator`](@ref), it is the *differential* entropy that is computed.

```@docs
EntropyEstimator
```

### Nearest neighbor based estimators

```@docs
Kraskov
KozachenkoLeonenko
Zhu
ZhuSingh
GaoNaive
GaoNaiveCorrected
Goria
LeonenkoProzantoSavani
```

### Order statistics based estimators

```@docs
Vasicek
AlizadehArghami
Ebrahimi
Correa
```

## Convenience

In this subsection we expand documentation strings of "entropy names" that are used commonly in the literature, such as "permutation entropy". As we made clear in [API & terminology](@ref), these are just the existing Shannon entropy with a particularly chosen probability estimator. We have only defined convenience functions for the most used names, and arbitrary more specialized convenience functions can be easily defined in a couple lines of code.

```@docs
entropy_permutation
entropy_spatial_permutation
entropy_wavelet
entropy_dispersion
```

## Examples

Here, we'll test the different nearest-neighbor based differential entropy estimators on a three-dimensional normal distribution
$\mathcal{N} (\mu, \Sigma)$ with zero means and covariance matrix $\Sigma = diag(r_1, r_2, r_3)$ with $r_1 = r_2 = r_3 = 0.5$. 
The analytical entropy for multivariate Gaussian is $H(\mathcal{N} (\mu, \Sigma)) = \dfrac{1}{2}\log(\det(2\pi e \Sigma))$. In our case, $\Sigma$ is diagonal, so $\det(\Sigma) = (0.5)^3$ and $H = 0.5\log(2\pi e (0.5)^3)\approx 3.217$.

Several of these estimators have been shown to convergence to the true entropy with an increasing number of samples. Therefore, we test the 
estimators on samples of increasing size $N$, where $N$ ranges from
1000 to 30000. Since we're estimating entropy from *samples* of
a normal distribution, we don't expect the estimates to perfectly match the analytical entropy every time.
On *average*, however, they should hit the target when the sample size
gets large enough.

### Analytical and estimated entropies

We'll first make two helper functions.

- **`analytical_entropy(estimators, Ls; d::Int, r, base = 2)`**: Computes the analytical  
Shannon differential entropy to the given `base` of a multivariate normal distribution with covariance matrix with diagonal elements `r` and zeros on the off-diagonal. Does so for each of the given `estimators` for each
sample size in `Ls`.
- **`mvnormal_entropies(; d::Int, r, base = 2, kwargs...)`**: Estimates  the Shannon entropy to the given `base` of samples from a multivariate normal distribution as specified as above.
- **`plot`

```@example ex_entropy_estimators
using CausalityTools
using Distributions: MvNormal
using LinearAlgebra
using Statistics: quantile
using Random; rng = MersenneTwister(12345678)
using CairoMakie

analytical_entropy(; d::Int, r, base = 2) = 
    0.5*log(det(2*pi*ℯ*diagm(repeat([r], d)))) / log(base, ℯ) # convert to desired base

function mvnormal_entropies(estimators, Ls; 
        d = 3,
        base = 2,
        nreps = 50,
        r = 0.5,
    )
    μ = zeros(d)
    Σ = diagm(repeat([r], d))
    N = MvNormal(μ, Σ)    
    Hs = [[zeros(nreps) for L in Ls] for est in estimators]
    data = [Dataset([rand(rng, N) for i = 1:maximum(Ls)]) for i = 1:nreps]
    for (e, est) in enumerate(estimators)
        for (l, L) in enumerate(Ls)
            for i = 1:nreps
                Hs[e][l][i] = entropy(Shannon(; base), est, data[i][1:L])
            end
        end
    end
    return Hs
end;
```

We'll also need a function to summarize the estimates.

```@example ex_entropy_estimators
# A helper to get the estimator name for plotting.
getname(est::EntropyEstimator) = typeof(est).name.name  |> string
function medians_and_quantiles(Hs, Ls; q = 0.95)
    medians = [zeros(length(Ls)) for est in estimators]
    lb = [zeros(length(Ls)) for est in estimators]
    ub = [zeros(length(Ls)) for est in estimators]

    for (e, est) in enumerate(estimators)
        for (l, L) in enumerate(Ls)
            ĥs = Hs[e][l] # nreps estimates for this combinations of e and l
            medians[e][l] = quantile(ĥs, 0.5)
            lb[e][l] = quantile(ĥs, (1 - q) / 2)
            ub[e][l] = quantile(ĥs, 1 - ((1 - q) / 2))
        end
    end

    return medians, lb, ub
end;
```

Now, make some plotting helper functions.

```@example ex_entropy_estimators
struct Cyclor{T} <: AbstractVector{T}
    c::Vector{T}
    n::Int
end
Cyclor(c) = Cyclor(c, 0)

Base.length(c::Cyclor) = length(c.c)
Base.size(c::Cyclor) = size(c.c)
Base.iterate(c::Cyclor, state=1) = Base.iterate(c.c, state)
Base.getindex(c::Cyclor, i) = c.c[(i-1)%length(c.c) + 1]
Base.getindex(c::Cyclor, i::AbstractArray) = c.c[i]
function Base.getindex(c::Cyclor)
    c.n += 1
    c[c.n]
end
Base.iterate(c::Cyclor, i = 1) = iterate(c.c, i)

COLORSCHEME = [
    "#D43F3AFF", "#EEA236FF", "#5CB85CFF", "#46B8DAFF",
    "#357EBDFF", "#9632B8FF", "#B8B8B8FF",
]

COLORS = Cyclor(COLORSCHEME)
LINESTYLES = Cyclor(string.(["--", ".-", ".", "--.", "---..."]))
MARKERS = Cyclor(string.([:circle, :rect, :utriangle, :dtriangle, :diamond,
    :pentagon, :cross, :xcross]))

function plot_entropy_estimates(Hs, Ls, Htrue)
    # Summarize data (medians[e][l]) is the median of the e-th estimator for the 
    # l-th sample size).
    medians, lbs, ubs = medians_and_quantiles(Hs, Ls);

    fig = Figure(resolution = (800, 1000))
    ymax = (vcat(Hs...) |> Iterators.flatten |> maximum) * 1.1
    ymin = (vcat(Hs...) |> Iterators.flatten |> minimum) * 0.9

    # We have 9 estimators, so place them on a 5-by-2 grid
    positions = (Tuple(c) for c in CartesianIndices((5, 2)))
    for (i, (est, c)) in enumerate(zip(estimators, positions))
        ax = Axis(fig[first(c), last(c)],
            xlabel = "Sample size (L)",
            ylabel = "Ĥ (bits)",
            title = getname(est)
        )
        ylims!(ax, (ymin, ymax))
        # Ground truth
        hlines!(ax, [Htrue], 
            linestyle = :dash, 
            color = :black,
            linewidth = 2,
        )
        # Estimates
        band!(ax, Ls, lbs[i], ubs[i], color = (COLORS[i], 0.5))
        lines!(ax, Ls, medians[i], 
            label = getname(est),
            linestyle = LINESTYLES[i],
            color = COLORS[i],
            marker = MARKERS[i],
            linewidth = 2
        )
    end
    fig
end;
```

Now, we can finally run an ensemble of tests and plot the
confidence bands against the ground truth. This

```@example ex_entropy_estimators
k = 4
estimators = [
    Kraskov(; k), 
    KozachenkoLeonenko(), 
    GaoNaive(; k),
    GaoNaiveCorrected(; k),
    ZhuSingh(; k),
    Zhu(; k),
    Goria(; k),
    LeonenkoProzantoSavani(; k),
    Lord(; k = k*5)
]
Ls = [100:100:1000 |> collect; 2500:2500:5000 |> collect]

d = 3
r = 0.5
nreps = 30
Hs = mvnormal_entropies(estimators, Ls; d, r, nreps)
Htrue = analytical_entropy(; d, r)
plot_entropy_estimates(Hs, Ls, Htrue)
```
