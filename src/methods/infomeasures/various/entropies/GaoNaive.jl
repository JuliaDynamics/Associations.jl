using Entropies: DifferentialEntropyEstimator
using StateSpaceSets: AbstractDataset, Dataset
using Neighborhood: Euclidean, KDTree, NeighborNumber, Theiler
using Neighborhood: bulksearch
using SpecialFunctions: digamma

export GaoNaive, GaoNaiveCorrected

"""
    GaoNaive <: DifferentialEntropyEstimator
    GaoNaive(k = 1, w = 0, base = 2)

The `GaoNaive` estimator (Gao et al., 2015) computes the [`Shannon`](@ref)
[`entropy`](@ref) to the given `base`, using a `k`-th nearest-neighbor approach
based on Singh (2003)[^Singh2003].

`w` is the Theiler window, which determines if temporal neighbors are excluded
during neighbor searches (defaults to `0`, meaning that only the point itself is excluded
when searching for neighbours).

[^Gao2005]:
    Gao, S., Ver Steeg, G., & Galstyan, A. (2015, February). Efficient estimation of
    mutual information for strongly dependent variables. In Artificial intelligence and
        statistics (pp. 277-286). PMLR.
[^Singh2003]:
    Singh, H., Misra, N., Hnizdo, V., Fedorowicz, A., & Demchuk, E. (2003). Nearest
    neighbor estimates of entropy. American journal of mathematical and management
    sciences, 23(3-4), 301-321.
"""
Base.@kwdef struct GaoNaive{B} <: DifferentialEntropyEstimator
    k::Int = 1
    w::Int = 0
    base::B = 2
end

function entropy(e::Renyi, est::GaoNaive, x::AbstractDataset{D}) where D
    e.q == 1 || throw(ArgumentError(
        "Renyi entropy with q = $(e.q) not implemented for $(typeof(est)) estimator"
    ))
    (; k, w) = est
    N = length(x)
    f = (k  * gamma(D / 2 + 1)) / ( (N - 1) * π^(D / 2))
    tree = KDTree(x, Euclidean())
    idxs, ds = bulksearch(tree, x, NeighborNumber(k), Theiler(w))
    h = -(1 / N) * sum(log(f * 1 / last(dᵢ)^D) for dᵢ in ds) # in nats
    return h / log(e.base, ℯ) # convert to target unit
end

"""
    GaoNaiveCorrected <: DifferentialEntropyEstimator
    GaoNaiveCorrected(k = 1, w = 0, base = 2)

The `GaoNaiveCorrected` estimator (Gao et al., 2015), computes the [`Shannon`](@ref)
[`entropy`](@ref) to the given `base`, using a `k`-th nearest-neighbor approach
based on Singh (2003)[^Singh2003].

`GaoNaiveCorrected` is identical to the [`GaoNaive`](@ref) estimator, except it adds
correction terms that ensures the estimator is asymptotically unbiased.

`w` is the Theiler window, which determines if temporal neighbors are excluded
during neighbor searches (defaults to `0`, meaning that only the point itself is excluded
when searching for neighbours).

[^Gao2005]:
    Gao, S., Ver Steeg, G., & Galstyan, A. (2015, February). Efficient estimation of
    mutual information for strongly dependent variables. In Artificial intelligence and
        statistics (pp. 277-286). PMLR.
[^Singh2003]:
    Singh, H., Misra, N., Hnizdo, V., Fedorowicz, A., & Demchuk, E. (2003). Nearest
    neighbor estimates of entropy. American journal of mathematical and management
    sciences, 23(3-4), 301-321.
"""
Base.@kwdef struct GaoNaiveCorrected{B} <: DifferentialEntropyEstimator
    k::Int = 1
    w::Int = 0
    base::B = 2
end

function entropy(e::Renyi, est::GaoNaiveCorrected, x::AbstractDataset{D}) where D
    (; k, w) = est
    N = length(x)
    f = (k  * gamma(D / 2 + 1)) / ( (N - 1) * π^(D / 2))
    tree = KDTree(x, Euclidean())
    idxs, ds = bulksearch(tree, x, NeighborNumber(k), Theiler(w))
    corr = digamma(k) - log(k)
    h = -(1 / N) * sum(log(f * 1 / last(dᵢ)^D) for dᵢ in ds) - corr # in nats
    return h / log(e.base, ℯ) # convert to target unit
end
