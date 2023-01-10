
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
