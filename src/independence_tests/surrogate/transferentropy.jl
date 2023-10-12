
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

function independence(test::SurrogateAssociationTest{<:TransferEntropy{<:E, <:EmbeddingTypes}}, x::AbstractVector...) where {E}
    (; measure, est, rng, surrogate, nshuffles) = test

    cmi = te_to_cmi(measure)
    Ŝ, T⁺, S, T, C = marginals_and_surrogenerator(measure.embedding, surrogate, x...; rng)
    TC = StateSpaceSet(T, C)
    @assert length(T⁺) == length(S) == length(TC)
    Î = estimate(cmi, est, T⁺, S, TC)
    Îs = zeros(nshuffles)
    for b in 1:nshuffles
        # TE(ŝ -> t) := I(t⁺; ŝ⁻ | t⁻, c⁻)
        Îs[b] = estimate(cmi, est, T⁺, Ŝ(), TC)
    end
    p = count(Î .<= Îs) / nshuffles

    return SurrogateAssociationTestResult(length(x), Î, Îs, p, nshuffles)
end

function independence(test::SurrogateAssociationTest{<:TransferEntropy{<:E, <:EmbeddingTypes}, <:TransferEntropyEstimator}, x::AbstractVector...) where {E}
    (; measure, est, rng, surrogate, nshuffles) = test

    Ŝ, T⁺, S, T, C = marginals_and_surrogenerator(measure.embedding, surrogate, x...; rng)
    @assert length(T⁺) == length(S) == length(T) == length(C)
    Î = estimate(measure, est, S, T, T⁺, C)
    Îs = zeros(nshuffles)
    for b in 1:nshuffles
        # TE(ŝ -> t) := I(t⁺; ŝ⁻ | t⁻, c⁻)
        Îs[b] = estimate(measure, est, StateSpaceSet(Ŝ()), T, T⁺, C)
    end
    p = count(Î .<= Îs) / nshuffles

    return SurrogateAssociationTestResult(length(x), Î, Îs, p, nshuffles)
end

function SurrogateAssociationTest(measure::TEShannon, est::Nothing, args...; kwargs...)
    txt = "A valid estimator must be provided as second argument to `SurrogateAssociationTest` " *
        "when using the `TEShannon` measure.\n" *
        "Do e.g. SurrogateAssociationTest(TEShannon(), FPVP())"
    throw(ArgumentError(txt))
end
