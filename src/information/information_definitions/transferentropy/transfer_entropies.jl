export TransferEntropy

"""
    TransferEntropy <: AssociationMeasure

The supertype of all transfer entropy measures. Concrete subtypes are
- [`TEShannon`](@ref)
- [`TERenyiJizba`](@ref)
"""
abstract type TransferEntropy <: MultivariateInformationMeasure end

max_inputs_vars(::TransferEntropy) = 3
is_directed(m::TransferEntropy) = true

include("embedding.jl")
include("utils/utils.jl")
include("utils/OptimiseTraditional.jl")

# If the estimator is not a dedicated `TransferEntropyEstimator`, then we need to 
# convert the estimator to a conditional mutual information estimator.
function information(est::MultivariateInformationMeasureEstimator{<:TransferEntropy}, x...)
    definition = est.definition
    # If a conditional input (x[3]) is not provided, then C is just a 0-dimensional
    # StateSpaceSet. The horizontal concatenation of C with T then just returns T.
    # We therefore don't need separate methods for the conditional and non-conditional
    # cases.
    S, T, T⁺, C = individual_marginals_te(definition.embedding, x...)

    cmi_est = convert_to_cmi_estimator(est)

    # Estimate by letting TE(s -> t | c) := I(t⁺; s⁻ | t⁻, c⁻). 
    return information(cmi_est, T⁺, S, StateSpaceSet(T, C))
end


function individual_marginals_te(emb::EmbeddingTE, x::AbstractVector...)
    joint, vars, τs, js = te_embed(emb, x...)
    S = joint[:, vars.S]
    T = joint[:, vars.T]
    Tf = joint[:, vars.Tf]
    C = joint[:, vars.C]
    return S, T, Tf, C
end

function h4_marginals(definition::TransferEntropy, x...)
    S, T, T⁺, C = individual_marginals_te(definition.embedding, x...)
    joint = StateSpaceSet(S, T, T⁺, C)
    ST = StateSpaceSet(S, T, C)
    TT⁺ = StateSpaceSet(T, T⁺, C)
    T = StateSpaceSet(T, C)
    return joint, ST, TT⁺, T
end

# Concrete implementations
include("TEShannon.jl")
include("TERenyiJizba.jl")


# Special estimation
include("transferoperator.jl")