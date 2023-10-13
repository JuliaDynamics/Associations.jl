export TERenyiJizba

"""
    TERenyiJizba() <: TransferEntropy

The Rényi transfer entropy from [Jizba2012](@citet).

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise
    and conditional dependence.
- Use with [`transferentropy`](@ref) to compute the raw transfer entropy.

## Description

The transfer entropy from source ``S`` to target ``T``, potentially
conditioned on ``C`` is defined as

```math
\\begin{align*}
TE(S \\to T) &:= I_q^{R_J}(T^+; S^- | T^-) \\\\
TE(S \\to T | C) &:= I_q^{R_J}(T^+; S^- | T^-, C^-),
\\end{align*},
```
where ``I_q^{R_J}(T^+; S^- | T^-)`` is Jizba et al. (2012)'s definition of
conditional mutual information ([`CMIRenyiJizba`](@ref)).
The variables ``T^+``, ``T^-``,
``S^-`` and ``C^-`` are described in the docstring for [`transferentropy`](@ref).

## Estimation

Estimating Jizba's Rényi transfer entropy is a bit complicated, since it doesn't have 
a dedicated estimator. Instead, we re-write the Rényi transfer entropy as a 
Rényi conditional mutual information, and estimate it using an 
[`EntropyDecomposition`](@ref) with a suitable discrete/differential Rényi entropy
estimator from the list below as its input.

| Estimator                      | Sub-estimator                    | Principle                    |
| :----------------------------- | :------------------------------- | :--------------------------- |
| [`EntropyDecomposition`](@ref) | [`LeonenkoProzantoSavani`](@ref) | Four-entropies decomposition |
| [`EntropyDecomposition`](@ref) | [`ValueBinning`](@ref)           | Four-entropies decomposition |
| [`EntropyDecomposition`](@ref) | [`Dispersion`](@ref)             | Four-entropies decomposition |
| [`EntropyDecomposition`](@ref) | [`OrdinalPatterns`](@ref)        | Four-entropies decomposition |
| [`EntropyDecomposition`](@ref) | [`UniqueElements`](@ref)         | Four-entropies decomposition |

Any of these estimators must be given as input to a [`CMIDecomposition](@ref) estimator.

## Example

```julia
using CausalityTools, Random
rng = MersenneTwister(1234)
x, y, z = rand(rng, 1000), rand(rng, 1000), rand(rng, 1000)
def = CMIRenyiJizba(q = 1.5)

# Using a differential Rényi entropy estimator
est = EntropyDecomposition(def, LeonenkoProzantoSavani(Renyi(), k = 10))
information(est, x, y, z)

# Using a plug-in Rényi entropy estimator, discretizing using ordinal patterns.
est = EntropyDecomposition(def, PlugIn(Renyi()), OrdinalPatterns(m=2), RelativeAmount())
information(est, x, y, z)
```
"""
struct TERenyiJizba{B, Q, EMB} <: TransferEntropy
    base::B
    q::Q
    embedding::EMB
    function TERenyiJizba(; base::B = 2, q::Q = 1.5, embedding::EMB = EmbeddingTE()) where {B, Q, EMB}
        return new{B, Q, EMB}(base, q, embedding)
    end
end

function convert_to_cmi_estimator(est::EntropyDecomposition{<:TERenyiJizba, <:DiscreteInfoEstimator{Renyi}})
    (; definition, est, discretization, pest) = est
    base = definition.base
    return EntropyDecomposition(CMIRenyiJizba(; base), est, discretization, pest)
end

function convert_to_cmi_estimator(est::EntropyDecomposition{<:TERenyiJizba, <:DifferentialInfoEstimator{Renyi}})
    return EntropyDecomposition(CMIRenyiJizba(; est.definition.base), est.est)
end

