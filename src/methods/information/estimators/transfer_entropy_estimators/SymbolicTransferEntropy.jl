export SymbolicTransferEntropy

# TODO: update to new syntax
"""
    SymbolicTransferEntropy <: TransferEntropyEstimator
    SymbolicTransferEntropy(definition = TEShannon(); m = 3, τ = 1, 
        lt = ComplexityMeasures.isless_rand

A convenience estimator for symbolic transfer entropy [Staniek2008](@cite).

## Compatible measures

- [`TEShannon`](@ref)

## Description

[Symbolic transfer entropy](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.158101)
consists of two simple steps. First, the input time series are encoded using [`codify`](@ref) with
the [`CodifyVariables`](@ref) discretization and the [`OrdinalPatterns`](@ref) outcome space. This 
transforms the input time series into integer time series. Transfer entropy entropy is then 
estimated from the encoded time series by applying  

Transfer entropy is then estimated as usual on the encoded timeseries with the embedding
dictated by `definition` and the [`JointProbabilities`](@ref) estimator.

## Examples

- [Example 1](@ref example_TEShannon_SymbolicTransferEntropy)
"""
struct SymbolicTransferEntropy{M} <: TransferEntropyEstimator{M}
    definition::M
    m::Int
    τ::Int
    lt::Function
end

function SymbolicTransferEntropy(definition::M = TEShannon(); 
        m = 3, τ = 1, lt = ComplexityMeasures.isless_rand) where M
    return SymbolicTransferEntropy{M}(definition, m, τ, lt)
end

function association(est::SymbolicTransferEntropy{<:TEShannon}, x::AbstractVector...)
    (; m, τ, lt) = est
    discretization = CodifyVariables(OrdinalPatterns(; m, τ, lt))

    x̂ = (codify(discretization, xᵢ) for xᵢ in x) 

    te_definition = est.definition
    embedding = te_definition.embedding
    # If a conditional input (x[3]) is not provided, then C is just a 0-dimensional
    # StateSpaceSet. The horizontal concatenation of C with T then just returns T.
    # We therefore don't need separate methods for the conditional and non-conditional
    # cases.
    S, T, T⁺, C = individual_marginals_te(embedding, x̂...)

    # We have already encoded the marginals, so when computing CMI, we can 
    # simply use `UniqueElements`.
    cmi_def = CMIShannon(; base = est.definition.base)
    disc = CodifyVariables(UniqueElements())
    
    est_unique = JointProbabilities(cmi_def, disc)
    return association(est_unique, T⁺, S,  StateSpaceSet(T, C))
end
