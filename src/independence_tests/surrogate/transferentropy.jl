
using TimeseriesSurrogates: Surrogate
using TimeseriesSurrogates: surrogenerator

const EmbeddingTypes = Union{EmbeddingTE, OptimiseTraditional}
function marginals_and_surrogenerator(emb::EmbeddingTE, surrogate::Surrogate, x::AbstractVector...; rng)
    if emb.dS > 1 && surrogate ∉ [RandomShuffle]
        throw(ArgumentError("Dimension of source embedding must be 1 to be applicable with surrogate methods"))
    end
    S, T, T⁺, C = individual_marginals_te(emb, x...)
    Ŝ = surrogenerator(S[:, 1], surrogate, rng)

    return Ŝ, T⁺, S, T, C
end
function marginals_and_surrogenerator(opt::OptimiseTraditional, surrogate::Surrogate, x::AbstractVector...; rng)
    emb = optimize_marginals_te(opt, x...; exclude_source = true)
    S, T, T⁺, C = individual_marginals_te(emb, x...)
    Ŝ = surrogenerator(S[:, 1], surrogate, rng)

    return Ŝ, T⁺, S, T, C
end

function independence(test::SurrogateAssociationTest{<:EntropyDecomposition{<:TransferEntropy}}, x, args...)
    (; est_or_measure, rng, surrogate, nshuffles) = test
    embedding = est_or_measure.definition.embedding

    cmi_est = convert_to_cmi_estimator(est_or_measure)
    Ŝ, T⁺, S, T, C = marginals_and_surrogenerator(embedding, surrogate, x, args...; rng)
    TC = StateSpaceSet(T, C)
    @assert length(T⁺) == length(S) == length(TC)
    Î = association(cmi_est, T⁺, S, TC)
    Îs = zeros(nshuffles)
    for b in 1:nshuffles
        # TE(ŝ -> t) := I(t⁺; ŝ⁻ | t⁻, c⁻)
        Îs[b] = association(cmi_est, T⁺, Ŝ(), TC)
    end
    p = count(Î .<= Îs) / nshuffles

    return SurrogateAssociationTestResult(length(x), Î, Îs, p, nshuffles)
end