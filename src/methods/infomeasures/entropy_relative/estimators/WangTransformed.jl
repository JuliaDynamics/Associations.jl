export WangTransformed

"""
    WangTransformed <: RelativeEntropyEstimator
    WangTransformedn(k = 5, l = 5, w = 0)

The `WangTransformed` relative entropy estimator (Wang et al., 2009[^Wang2009] is
identical to [`Wang`](@ref), but proprocesses the input data so that
their covariance matrix is close to the identity matrix.

[^Wang2009]:
    Wang, Q., Kulkarni, S. R., & Verdú, S. (2009). Divergence estimation for
    multidimensional densities via k-Nearest-Neighbor distances. IEEE Transactions on
    Information Theory, 55(5), 2392-2405.
"""
Base.@kwdef struct WangTransformed <: RelativeEntropyEstimator
    k::Int = 5
    l::Int = 5
    w::Int = 0
end

function transform_samples(est::WangTransformed, x, y)
    D = append!(copy(x), y)
    μ̂ = 1/length(D) * mean(D)
    # TODO: explicit loop will be faster because we don't need to compute differences twice
    Ĉ = 1/(length(D)-1) * sum((x̂ᵢ - μ̂)*transpose(x̂ᵢ - μ̂) for x̂ᵢ in D)
    Ĉinvsq = 1 ./ sqrt(Ĉ)
    D̂ = Dataset([Ĉinvsq * (dᵢ - μ̂) for dᵢ in D])
end

function entropy_relative(e::Renyi, est::WangTransformed,
        x::AbstractDataset{D}, y::AbstractDataset{D}) where D
    (; k, l, w) = est
    D̂ = transform_samples(est, x, y)
    N = length(D̂) ÷ 2
    X̂ = D̂[1:N]
    Ŷ = D̂[(N + 1):end]
    est_untransformed = Wang(; k, l, w)

    return entropy_relative(e, est_untransformed, X̂, Ŷ)
end

entropy_relative(est::WangTransformed, args...; base = 2, kwargs...) =
    entropy(Shannon(; base), est, args...; kwargs...)
