# [Entropy (generalized)](@id entropies)

## API & entropy types

```@docs
entropy
Entropy
```

```@docs
Shannon
Renyi
Tsallis
Curado
StretchedExponential
```

## Estimators

Discrete entropies are just simple functions (sums, actually) of 
probability mass functions [(pmf)](https://en.wikipedia.org/wiki/Probability_mass_function),
which you can estimate using [`ProbabilitiesEstimator`](@ref)s. Every 
[`ProbabilitiesEstimator`](@ref)s thus serves as a naive plug-in estimators of all 
generalized entropies. No bias correction is currently applied to any of the discrete 
estimators.

Continuous (differential) entropies, on the other hand, are approximated by [`EntropyEstimator`](@ref)s,
and are functions of *integrals*, and rely on estimating density functionals. 
Each [`EntropyEstimator`](@ref)s uses a different technique for approximating relevant
probability densities, and is often tailed to one or a few types of generalized entropy.
For example, [`Kraskov`](@ref) estimates the [`Shannon`](@ref) entropy, while
[`LeonenkoProzantoSavani`](@ref) estimates [`Shannon`](@ref), [`Renyi`](@ref), and 
[`Tsallis`](@ref) entropies.

```@docs
EntropyEstimator
```

```@docs
Kraskov
KozachenkoLeonenko
Zhu
ZhuSingh
GaoNaive
GaoNaiveCorrected
Goria
LeonenkoProzantoSavani
Vasicek
AlizadehArghami
Ebrahimi
Correa
```

## Implementations

The table below summarizes this available implementations. An `x` indicates that no
implementation exists. The table scrolls sideways, so is best viewed on a large screen.

| Estimator                                   | Principle                   | Input data          | [`Shannon`](@ref) | [`Renyi`](@ref) | [`Tsallis`](@ref) | [`Curado`](@ref) | [`StretchedExponential`](@ref) |
| ------------------------------------------- | --------------------------- | ------------------- | :---------------: | :-------------: | :---------------: | :--------------: | :----------------------------: |
| [`CountOccurrences`](@ref)                  | Frequencies                 | `Vector`, `Dataset` |     Discrete      |    Discrete     |     Discrete      |     Discrete     |            Discrete            |
| [`ValueHistogram`](@ref)                    | Binning (histogram)         | `Vector`, `Dataset` |     Discrete      |    Discrete     |     Discrete      |     Discrete     |            Discrete            |
| [`TransferOperator`](@ref)                  | Binning (transfer operator) | `Vector`, `Dataset` |     Discrete      |    Discrete     |     Discrete      |     Discrete     |            Discrete            |
| [`NaiveKernel`](@ref)                       | Kernel density estimation   | `Dataset`           |     Discrete      |    Discrete     |     Discrete      |     Discrete     |            Discrete            |
| [`LocalLikelihood`](@ref)                   | Local likelhood Estimation  | `Dataset`           |     Discrete      |    Discrete     |     Discrete      |     Discrete     |            Discrete            |
| [`SymbolicPermutation`](@ref)               | Ordinal patterns            | `Vector`            |     Discrete      |    Discrete     |     Discrete      |     Discrete     |            Discrete            |
| [`SymbolicWeightedPermutation`](@ref)       | Ordinal patterns            | `Vector`            |     Discrete      |    Discrete     |     Discrete      |     Discrete     |            Discrete            |
| [`SymbolicAmplitudeAwarePermutation`](@ref) | Ordinal patterns            | `Vector`            |     Discrete      |    Discrete     |     Discrete      |     Discrete     |            Discrete            |
| [`Dispersion`](@ref)                        | Dispersion patterns         | `Vector`            |     Discrete      |    Discrete     |     Discrete      |     Discrete     |            Discrete            |
| [`Diversity`](@ref)                         | Cosine similarity           | `Vector`            |     Discrete      |    Discrete     |     Discrete      |     Discrete     |            Discrete            |
| [`WaveletOverlap`](@ref)                    | Wavelet transform           | `Vector`            |     Discrete      |    Discrete     |     Discrete      |     Discrete     |            Discrete            |
| [`PowerSpectrum`](@ref)                     | Fourier spectra             | `Vector`, `Dataset` |     Discrete      |    Discrete     |     Discrete      |     Discrete     |            Discrete            |
| [`KozachenkoLeonenko`](@ref)                | Nearest neighbors           | `Dataset`           |    Continuous     |        x        |         x         |        x         |               x                |
| [`Kraskov`](@ref)                           | Nearest neighbors           | `Dataset`           |    Continuous     |        x        |         x         |        x         |               x                |
| [`Zhu`](@ref)                               | Nearest neighbors           | `Dataset`           |    Continuous     |        x        |         x         |        x         |               x                |
| [`ZhuSingh`](@ref)                          | Nearest neighbors           | `Dataset`           |    Continuous     |        x        |         x         |        x         |               x                |
| [`GaoNaive`](@ref)                          | Nearest neighbors           | `Dataset`           |    Continuous     |        x        |         x         |        x         |               x                |
| [`GaoNaiveCorrected`](@ref)                 | Nearest neighbors           | `Dataset`           |    Continuous     |        x        |         x         |        x         |               x                |
| [`Goria`](@ref)                             | Nearest neighbors           | `Dataset`           |    Continuous     |        x        |         x         |        x         |               x                |
| [`LeonenkoProzantoSavani`](@ref)            | Nearest neighbors           | `Dataset`           |    Continuous     |   Continuous    |    Continuous     |        x         |               x                |
| [`Vasicek`](@ref)                           | Order statistics            | `Vector`            |    Continuous     |        x        |         x         |        x         |               x                |
| [`Ebrahimi`](@ref)                          | Order statistics            | `Vector`            |    Continuous     |        x        |         x         |        x         |               x                |
| [`Correa`](@ref)                            | Order statistics            | `Vector`            |    Continuous     |        x        |         x         |        x         |               x                |
| [`AlizadehArghami`](@ref)                   | Order statistics            | `Vector`            |    Continuous     |        x        |         x         |        x         |               x                |

!!! info "Missing implementation?"
    Submit an issue or a pull request! We're happy to assist.

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
