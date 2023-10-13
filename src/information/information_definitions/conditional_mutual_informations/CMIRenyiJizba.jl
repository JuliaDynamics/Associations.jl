using ComplexityMeasures: Renyi

export CMIRenyiJizba

"""
    CMIRenyiJizba <: ConditionalMutualInformation
    CMIRenyiJizba(; base = 2, q = 1.5)

The Rényi conditional mutual information ``I_q^{R_{J}}(X; Y | Z`` defined in
[Jizba2012](@citet).

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise dependence.
- Use with [`condmutualinfo`](@ref) to compute the raw conditional mutual information.

## Definition

```math
I_q^{R_{J}}(X; Y | Z) = I_q^{R_{J}}(X; Y, Z) - I_q^{R_{J}}(X; Z),
```

where ``I_q^{R_{J}}(X; Z)`` is the [`MIRenyiJizba`](@ref) mutual information.

## Estimation

The following estimators can be used to compute `CMIRenyiJizba`.

| Estimator                      |          Sub-estimator           | Principle                    |
| ------------------------------ | :------------------------------: | :--------------------------- |
| [`JointProbabilities`](@ref)   |                -                 | Discrete joint pmf           |
| [`EntropyDecomposition`](@ref) | [`LeonenkoProzantoSavani`](@ref) | Four-entropies decomposition |

## Examples

```julia
using CausalityTools, Random
rng = MersenneTwister(1234)
x, y, z = rand(rng, 1000), rand(rng, 1000), rand(rng, 1000)
def = CMIRenyiJizba(q = 1.5)
est = EntropyDecomposition(def, LeonenkoProzantoSavani(Renyi(), k = 10))
information(est, x, y, z)

est = EntropyDecomposition(def, PlugIn(Renyi()), ValueBinning(3), RelativeAmount())
information(est, x, y, z)
```
"""
Base.@kwdef struct CMIRenyiJizba{B, Q} <: ConditionalMutualInformation
    base::B = 2
    q::Q = 1.5
end

function information(est::JointProbabilities{<:CMIRenyiJizba}, x, y, z)
    throw(ArgumentError("CMIRenyiJizba not implemented for $(typeof(est).name.name)"))
end

# --------------------------------------------------------------
# Conditional mutual information through entropy decomposition
# --------------------------------------------------------------
function information(est::EntropyDecomposition{<:CMIRenyiJizba, <:DifferentialInfoEstimator{<:Renyi}}, x, y, z)
    HXZ, HYZ, HXYZ, HZ = marginal_entropies_cmi4h_differential(est, x, y, z)
    cmi = HXZ + HYZ - HXYZ - HZ
    return cmi
end

function information(est::EntropyDecomposition{<:CMIRenyiJizba, <:DiscreteInfoEstimator{<:Renyi}}, x, y, z)
    HXZ, HYZ, HXYZ, HZ = marginal_entropies_cmi4h_discrete(est, x, y, z)
    cmi = HXZ + HYZ - HXYZ - HZ
    return cmi
end

# ------------------------------------------------
# Pretty printing for decomposition estimators.
# ------------------------------------------------
function decomposition_string(
        definition::CMIRenyiJizba, 
        est::EntropyDecomposition{<:CMIRenyiJizba, <:DiscreteInfoEstimator{<:Renyi}}
    ) 
    return "CMI_RJ(X, Y | Z) = H_R(X,Z) + H_R(Y,Z) - H_R(X,Y,Z) - H_R(Z)";
end

function decomposition_string(
    definition::CMIRenyiJizba, 
    est::EntropyDecomposition{<:CMIRenyiJizba, <:DifferentialInfoEstimator{<:Renyi}}
) 
    return "CMI_RJ(X, Y | Z) = h_R(X,Z) + h_R(Y,Z) - h_R(X,Y,Z) - h_R(Z)";
end
