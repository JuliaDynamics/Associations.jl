using Neighborhood: Chebyshev, KDTree, Theiler, NeighborNumber
using Neighborhood: bulksearch
using Distances: evaluate
using DelayEmbeddings.StateSpaceSets: SubDataset
using LinearAlgebra: det

"""
    Gao2017 <: EntropyEstimator
    Gao2017(k = 1, w = 1)

The `Gao2017` estimator (Gao et al. ,2017)[^Gao2017] can be used to estimate
[`Renyi`](@ref) differential entropy.

It does so by considering, for ``\\alpha \\neq 1``  the integral

```math
J_{\\alpha}(X) = \\int_{\\mathbb{R}^D} f^{\\alpha}(x) dx
```


[^Gao2017]: Gao, W., Oh, S., & Viswanath, P. (2017, June). Density functional estimators
    with k-nearest neighbor bandwidths. In 2017 IEEE International Symposium on Information
    Theory (ISIT) (pp. 1351-1355). IEEE.
"""
Base.@kwdef struct Gao2017{B, M} <: InformationEstimator
    k::Int = 1
    w::Int = 0
    metric::M = Euclidean()
end

function Î(q, est::Gao2017, x::AbstractDataset{D}) where D
    (; k, w, metric) = est
    N = length(x)
    tree = KDTree(x, metric)
    Bk,d,α,K = bias(est)
    idxs, ds = bulksearch(tree, x, NeighborNumber(k), Theiler(w))

end

# TODO: implement
multiplicative_bias(est::Gao2017) = 1.0

Base.@kwdef struct LocalLikelihood{M} <: ProbabilitiesEstimator
    k::Int = 1
    w::Int = 0
    metric::M = Euclidean()
end

function point_densities(est::LocalLikelihood, x::AbstractDataset{D}) where D

    (; k, w, metric) = est

    N = length(x)
    kmax = max(floor(Int, log(N)), k)

    tree = KDTree(x)
    # The bandwidth `ds[i]` for the point `x[i]` is the distance to the `k`-th nearest
    # neighbor of `x[i]`.
    idxs, ds = bulksearch(tree, x, NeighborNumber(kmax), Theiler(w))
    hs = [d[k] for d in ds]

    densities = zeros(N)
    for i = 1:N
        xᵢ = x[i]
        hᵢ = hs[i]
        neighborsᵢ = @view x[idxs[i]]
        densities[i] = point_density(est, xᵢ, hᵢ, neighborsᵢ)
    end
    return densities
end

# Compute the local density around point xᵢ, given its `neighborsᵢ`
function point_density(est::LocalLikelihood, xᵢ, hᵢ, neighborsᵢ::SubDataset{D}) where D
    S₀ = 0.0
    S₁ = zeros(MVector{D, Float64})
    S₂ = zeros(MMatrix{D, D, Float64})
    # Gao et al, in the original paper, only loops through the floor(Int, log(N)) nearest
    # neighbors of x[i]. No need to go through all.
    hᵢsq = hᵢ^2
    twicehᵢsq = 2*hᵢsq
    for (k, nⱼ) in enumerate(neighborsᵢ)
        dᵢ = evaluate(est.metric, xᵢ, nⱼ)
        eᵢ = exp(-dᵢ / twicehᵢsq)
        xdiff = (nⱼ - xᵢ)
        S₀ += eᵢ
        S₁ .+= xdiff * (eᵢ / hᵢ)
        S₂ .+= (xdiff * transpose(xdiff)) .* (eᵢ / twicehᵢsq)
    end

    μ = S₁ / S₀
    Σ = (S₂ / S₀) - (S₁ * transpose(S₁) / (S₀^2))

    num = (S₀ * exp(-0.5*transpose(μ) * inv(Σ) * μ))
    den = (N * (2π)^(D/2) * (hᵢ^D) * det(Σ)^(1/2))
    return num / den
end
