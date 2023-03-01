using Neighborhood: Euclidean, Chebyshev, KDTree, Theiler, NeighborNumber
using Neighborhood: bulksearch
using Distances: evaluate
using DelayEmbeddings.StateSpaceSets: SubStateSpaceSet
using LinearAlgebra: det, norm
using StateSpaceSets: StateSpaceSet
using StaticArrays: MVector, MMatrix, SVector, SMatrix, @SVector

import ComplexityMeasures: entropy, total_outcomes, outcomes, probabilities, probabilities_and_outcomes

export LocalLikelihood
"""
    LocalLikelihood <: ProbabilitiesEstimator
    LocalLikelihood(k = 5, w = 0, metric = Euclidean())

The `LocalLikelihood` estimator estimates the density around a given query point
by a Gaussian kernel informed by the local mean and covariance.

To form probabilities from the pointwise density estimates, the densities are
simply sum-normalized to 1.

## Outcome space

The [`outcome_space`](@ref) for `LocalLikelihood` is the set of input points.
"""
Base.@kwdef struct LocalLikelihood{M} <: ProbabilitiesEstimator
    k::Int = 5
    w::Int = 0
    metric::M = Euclidean()
end

function point_densities(est::LocalLikelihood, x::AbstractStateSpaceSet{D}) where D
    (; k, w, metric) = est
    N = length(x)
    # Modified heuristic from Gao et al. (2017): it is sufficient to consider the
    # `K = max(floor(Int, log(N), k)` nearest neighbors neighbors of `x[i]` when
    # estimating the local density. A global point-search is pointless and expensive.
    kmax = max(floor(Int, log(N)), k)

    # The bandwidth `bws[i]` for the point `x[i]` is the distance to the `k`-th nearest
    # neighbor of `x[i]`. The local density around, in contrast, in formed by the `kmax`
    # nearest neighbors.
    tree = KDTree(x, Euclidean())
    idxs, ds = bulksearch(tree, x, NeighborNumber(kmax), Theiler(w))
    bws = [d[k] for d in ds]

    S₁ = zeros(MVector{D, Float64})
    S₂ = zeros(MMatrix{D, D, Float64})
    densities = zeros(N)
    for i = 1:N
        xᵢ = x[i]
        bwᵢ = bws[i]
        neighborsᵢ = x[idxs[i]]
        densities[i] = point_density!(S₁, S₂, est, xᵢ, bwᵢ, neighborsᵢ, N)
    end

    return densities
end

"""
    point_density!(S₁, S₂, est::LocalLikelihood, xᵢ, bwᵢ,
        neighborsᵢ::AbstractStateSpaceSet{D}) where D

Estimate the density around point `xᵢ` using a local likehood estimator, which is
a generalization of kernel density estimation. This is done by fitting a local gaussian
distribution around `xᵢ` from its local neighborhood (represented the points `neighborsᵢ`).
The bandwidth  `bwᵢ` is given by the distance from `xᵢ` to its `k`-th nearest neighbor.

`S₁` is a pre-allocated length-`D` vector which holds the means, and `S₂` is a pre-allocated
`D`-by-`D` matrix which holds the covariances. Both `S₁` and `S₂` are zeroed every time
`point_density!` is called.
"""
function point_density!(S₁, S₂, est::LocalLikelihood, xᵢ, bwᵢ, neighborsᵢ::AbstractStateSpaceSet{D}, Ntot::Int) where D
    N = length(neighborsᵢ)
    S₀ = 0.0;
    S₁ .= 0.0
    S₂ .= 0.0
    bwᵢ_sq = bwᵢ^2
    twice_bwᵢ_sq = 2*bwᵢ_sq
    for (k, nⱼ) in enumerate(neighborsᵢ)
        Δⱼ = (nⱼ - xᵢ)
        dᵢ = evaluate(est.metric, nⱼ, xᵢ)^2
        eᵢ = exp(-dᵢ / twice_bwᵢ_sq)
        S₀ += eᵢ
        S₁ += eᵢ * (Δⱼ / bwᵢ)
        S₂ += eᵢ * ((Δⱼ * transpose(Δⱼ)) / bwᵢ_sq)
    end
    # Weighted sample mean and sample variance
    μ = S₁ / S₀
    Σ = S₂ / S₀ - S₁*transpose(S₁) / S₀^2

    # if Σ is singular, we can't take its inverse either, so just return 0.0
    # density straight away. Heuristic from origina paper.
    detΣ = det(Σ)
    if detΣ < 1e-4^D
        return 0.0
    end

    # The commented-out code follows the paper. This gives nonsense results.
    #num = S₀ * exp(-0.5 * transpose(μ) * inv(Σ) * μ)
    #den = N * (2π)^(D/2) * (bwᵢ^D) * sqrt(detΣ)
    #return #num/den
    # the following code is from https://github.com/wgao9/lnn/blob/master/lnn.py,
    # by one of the original authors,
    # but I have no idea where this formula comes from. It seems to work, though.
    offset = transpose(μ) * inv(Σ) * μ
    return -log(S₀) +
        log(Ntot - 1) +
        0.5*D*log(2π) +
        D*log(bwᵢ) +
        0.5*log(detΣ) + 0.5*offset[1, 1]
end

function probabilities_and_outcomes(est::LocalLikelihood, x)
    return Probabilities(point_densities(est, x)), x
end
probabilities(est::LocalLikelihood, x) = Probabilities(point_densities(est, x))
outcomes(est::LocalLikelihood, x) = x
total_outcomes(x, est::LocalLikelihood) = length(x)

function entropy(e::Renyi, est::LocalLikelihood, x)
    !(e.q ≈ 1.0) || error("Renyi entropy for $(typeof(est)) estimator not defined for q = $(e.q) (i.e. Shannon entropy not defined)")
    N = length(x)
    ρs = point_densities(est, x)
    ĴkLNN = sum(ρs .^ (e.q - 1)) / (bias(e, est, x) * N)
    h = ĴkLNN / (e.q - 1)
    return _convert_logunit(h, ℯ, e.base)
end

function pt_in_unit_sphere(dim::Int)
    u = @SVector randn(dim)
    c = rand()^(1/dim)
    m = sqrt(sum(u .^ 2))
    v = u ./ m .* c
    return v
end
pts_in_unit_sphere(dim::Int, n::Int) = StateSpaceSet([pt_in_unit_sphere(dim) for i = 1:n])


# TODO: implement. not sure how, though. Gao (2017) is not very clear...
bias(e::Renyi, est::LocalLikelihood, x) = 1.0
