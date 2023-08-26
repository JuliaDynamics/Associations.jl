
Base.@deprecate te_embed(x, y, emb::EmbeddingTE) te_embed(emb::EmbeddingTE, x, y)
Base.@deprecate te_embed(x, y, z, emb::EmbeddingTE) te_embed(emb::EmbeddingTE, x, y, z)

function transferentropy(x::AbstractVector, y::AbstractVector, est; kwargs...)
    @warn """`transferentropy(x, y, est)` is deprecated. Use \
        `transferentropy(TEShannon(), est, x, y) instead`."""
    emb = EmbeddingTE(; kwargs...)
    return transferentropy(TEShannon(; embedding = emb, base = 2), est, x, y)
end

function transferentropy(x::AbstractVector, y::AbstractVector, z::AbstractVector,
        est; kwargs...)
    @warn """`transferentropy(x, y, est)` is deprecated. Use \
        `transferentropy(TEShannon(), est, x, y) instead`."""
    emb = EmbeddingTE(; kwargs...)
    return transferentropy(TEShannon(; embedding = emb, base = 2), est, x, y, z)
end

# Base.@deprecate bbnue(
#     source::AbstractVector, target::AbstractVector,
#     est::Union{ProbabilitiesEstimator, DifferentialInformationEstimator};
#     base = 2, η = 1,
#     include_instantaneous = true,
#     method_delay = "ac_min",
#     maxlag::Union{Int, Float64} = 0.05,
#     α = 0.05, nsurr = 19,
#     surr::TimeseriesSurrogates.Surrogate = RandomShuffle()
#     ) bbnue(TEShannon(; base), est, source, target; η, include_instantaneous,
#         method_delay, maxlag, α, nsurr, surr)
# Base.@deprecate bbnue(source::AbstractVector, target::AbstractVector, cond::AbstractVector,
#     est::Union{ProbabilitiesEstimator, DifferentialInformationEstimator};
#     η, include_instantaneous, method_delay, maxlag, α, nsurr,
#     surr) bbnue(TEShannon(; base), est, source, target, cond; η, include_instantaneous,
#         method_delay, maxlag, α, nsurr, surr)
