import Statistics.cov
using StaticArrays: SMatrix, @MMatrix, @MVector

"""
    Gao2017 <: DifferentialEntropyEstimator
    Gao2017(k = 1, w = 1, base = 2)

A resubstitution Shannon entropy estimator from [Gao2017](@citet).
"""
Base.@kwdef struct Gao2017{B, M} #<: CausalityTools.InformationEstimator
    k::Int = 1
    w::Int = 0
    metric::M = Euclidean()
end

function Î(q, est::Gao2017, x::AbstractStateSpaceSet{D}) where D
    (; k, w, metric) = est
    N = length(x)
    tree = KDTree(x, metric)
    Bk,d,α,K = bias(est)
    idxs, ds = bulksearch(tree, x, NeighborNumber(k), Theiler(w))

    # In the case of a multivariate Gaussian, maximum likelihood estimation simply
    # amounts to finding the sample means and sample covariance matrices.
    for xᵢ in x
        μ, Σ = mean_and_cov(x)
    end
end

import Statistics.cov

# Non-allocating and more than twice as fast as writing a wrapper
# `f(x) = Statistics.cov(Matrix(x))`.
# Also accepts SubStateSpaceSets, so we can use views on neighbor points.
function cov(x̄, x::AbstractStateSpaceSet{D}) where D
    N = length(x) - 1
    C = @MMatrix zeros(D, D)
    x̄ = mean(x)
    Δx = @MVector zeros(D)
    @inbounds for xᵢ in x
        Δx .= xᵢ - x̄
        C .+= Δx * transpose(Δx)
    end
    C ./= N
    return SMatrix{D, D}(C)
end
# So we don't have to compute the mean twice at every iteration.
cov(x::AbstractStateSpaceSet{D}) where D = cov(mean(x), x)
function mean_and_cov(x::AbstractStateSpaceSet{D}) where D
    μ = mean(x)
    Σ = cov(μ, x)
    return μ, Σ
end

# TODO: implement
multiplicative_bias(est::Gao2017) = 1.0
