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


# If the estimator is not a dedicated `TransferEntropyEstimator`, then we
# convert the estimator to a conditional mutual information estimator which we apply
# to appropriately constructed marginals.
function association(est::MultivariateInformationMeasureEstimator{<:TransferEntropy}, x...)
    te_definition = est.definition
    embedding = te_definition.embedding
    # If a conditional input (x[3]) is not provided, then C is just a 0-dimensional
    # StateSpaceSet. The horizontal concatenation of C with T then just returns T.
    # We therefore don't need separate methods for the conditional and non-conditional
    # cases.
    S, T, T⁺, C = individual_marginals_te(embedding, x...)

    cmi_est = convert_to_cmi_estimator(est)
    if est isa JointProbabilities
        # If discrete, we must codify each marginal separately and 
        # collect the "marginal marginals" into new statespace sets.
        tmp_TC = codify(est.discretization, StateSpaceSet(T, C))
        if tmp_TC isa AbstractVector
            T̂Ĉ = StateSpaceSet(tmp_TC)
        else
            T̂Ĉ = StateSpaceSet(tmp_TC...,)
        end
        Ŝ = codify(est.discretization, S)
        T̂⁺ = codify(est.discretization, T⁺)
        # We have already encoded the marginals, so when computing CMI, we can 
        # simply use `UniqueElements`.
        disc = CodifyVariables(UniqueElements())
        est_unique = JointProbabilities(cmi_est.definition, disc, est.pest)
        return association(est_unique, T̂⁺, Ŝ , T̂Ĉ)

    else
        #Estimate by letting TE(s -> t | c) := I(t⁺; s⁻ | t⁻, c⁻). 
        return association(cmi_est, T⁺, S, StateSpaceSet(T, C))
    end

end

function individual_marginals_te(emb::EmbeddingTE, x::VectorOr1DDataset...)
    joint, vars, τs, js = te_embed(emb, x...)
    S = joint[:, vars.S]
    T = joint[:, vars.T]
    Tf = joint[:, vars.Tf]
    C = joint[:, vars.C]
    return S, T, Tf, C
end

# function h4_marginals(definition::TransferEntropy, x...)
#     S, T, T⁺, C = individual_marginals_te(definition.embedding, x...)
#     joint = StateSpaceSet(S, T, T⁺, C)
#     ST = StateSpaceSet(S, T, C)
#     TT⁺ = StateSpaceSet(T, T⁺, C)
#     T = StateSpaceSet(T, C)
#     return joint, ST, TT⁺, T
# end

# Concrete implementations
include("TEShannon.jl")
include("TERenyiJizba.jl")


# Special estimation
include("transferoperator.jl")