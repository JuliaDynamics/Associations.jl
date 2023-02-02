
using TimeseriesSurrogates: Surrogate
using TimeseriesSurrogates: surrogenerator

const EmbeddingTypes = Union{EmbeddingTE, OptimiseTraditional}
function marginals_and_surrogenerator(emb::EmbeddingTE, surrogate::Surrogate, x::AbstractVector...; rng)
    if emb.dS != 1
        throw(ArgumentError("Dimension of source embedding must be 1 to be applicable with surrogate methods"))
    end
    S, T, T⁺, C = individual_marginals_te(emb, x...)
    Ŝ = surrogenerator(S[:, 1], surrogate, rng)

    return Ŝ, T⁺, S, Dataset(T, C)
end
function marginals_and_surrogenerator(opt::OptimiseTraditional, surrogate::Surrogate, x::AbstractVector...; rng)
    emb = optimize_marginals_te(opt, x...; exclude_source = true)
    S, T, T⁺, C = individual_marginals_te(emb, x...)
    Ŝ = surrogenerator(S[:, 1], surrogate, rng)

    return Ŝ, T⁺, S, Dataset(T, C)
end

function independence(test::SurrogateTest{<:TransferEntropy{<:E, <:EmbeddingTypes}}, x::AbstractVector...) where {E}
    (; measure, est, rng, surrogate, nshuffles) = test

    cmi = te_to_cmi(measure)
    Ŝ, T⁺, S, T = marginals_and_surrogenerator(measure.embedding, surrogate, x...; rng)
    @assert length(T⁺) == length(S) == length(T)
    Î = estimate(cmi, est, T⁺, S, T)
    Îs = zeros(nshuffles)
    for b in 1:nshuffles
        # TE(ŝ -> t) := I(t⁺; ŝ⁻ | t⁻, c⁻).
        Îs[b] = estimate(cmi, est, T⁺, Ŝ(), T)
    end
    p = count(Î .<= Îs) / nshuffles

    return SurrogateTestResult(Î, Îs, p, nshuffles)
end
