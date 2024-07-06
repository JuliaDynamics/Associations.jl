using ComplexityMeasures: Renyi

export CMIRenyiJizba

"""
    CMIRenyiJizba <: ConditionalMutualInformation
    CMIRenyiJizba(; base = 2, q = 1.5)

The Rényi conditional mutual information ``I_q^{R_{J}}(X; Y | Z)`` defined in
[Jizba2012](@citet).

## Usage

- Use with [`association`](@ref) to compute the raw  Rényi-Jizba conditional mutual information
    using of of the estimators listed below.
- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise conditional 
    independence using the Rényi-Jizba conditional mutual information.

## Compatible estimators

- [`JointProbabilities`](@ref)
- [`EntropyDecomposition`](@ref)

## Definition

```math
I_q^{R_{J}}(X; Y | Z) = I_q^{R_{J}}(X; Y, Z) - I_q^{R_{J}}(X; Z),
```

where ``I_q^{R_{J}}(X; Z)`` is the [`MIRenyiJizba`](@ref) mutual information.

## Examples

```julia
using CausalityTools, Random
rng = MersenneTwister(1234)
x, y, z = rand(rng, 1000), rand(rng, 1000), rand(rng, 1000)
def = CMIRenyiJizba(q = 1.5)

# Using a differential Rényi entropy estimator
est = EntropyDecomposition(def, LeonenkoProzantoSavani(Renyi(), k = 10))
association(est, x, y, z)

# Using a plug-in Rényi entropy estimator, discretizing using ordinal patterns.
est = EntropyDecomposition(def, PlugIn(Renyi()), OrdinalPatterns(m=2), RelativeAmount())
association(est, x, y, z)
```
"""
Base.@kwdef struct CMIRenyiJizba{B, Q} <: ConditionalMutualInformation
    base::B = 2
    q::Q = 1.5
end

function association(est::JointProbabilities{<:CMIRenyiJizba}, x, y, z)
    throw(ArgumentError("CMIRenyiJizba not implemented for $(typeof(est).name.name)"))
end

# --------------------------------------------------------------
# Conditional mutual information through entropy decomposition
# --------------------------------------------------------------
function association(est::EntropyDecomposition{<:CMIRenyiJizba, <:DifferentialInfoEstimator{<:Renyi}}, x, y, z)
    HXZ, HYZ, HXYZ, HZ = marginal_entropies_cmi4h_differential(est, x, y, z)
    cmi = HXZ + HYZ - HXYZ - HZ
    return cmi
end

function association(est::EntropyDecomposition{<:CMIRenyiJizba, <:DiscreteInfoEstimator{<:Renyi}}, x, y, z)
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
    return "Iᵣⱼ(X, Y | Z) = Hᵣ(X,Z) + Hᵣ(Y,Z) - Hᵣ(X,Y,Z) - Hᵣ(Z)";
end

function decomposition_string(
    definition::CMIRenyiJizba, 
    est::EntropyDecomposition{<:CMIRenyiJizba, <:DifferentialInfoEstimator{<:Renyi}}
) 
    return "Iᵣⱼ(X, Y | Z) = hᵣ(X,Z) + hᵣ(Y,Z) - hᵣ(X,Y,Z) - hᵣ(Z)";
end

# ---------------------------------
# Avoid some common errors
# ---------------------------------
function verify_decomposition_entropy_type(definition::CMIRenyiJizba, est::ENTROPY_ESTS)
    if !(est.definition isa Renyi)
        T = typeof(est.definition).name.name
        msg = "Can't decompose CMIRenyiJizba into a combination of $T entropies. Please provide a `Renyi` entropy estimator instead."
        throw(ArgumentError(msg))
    end
end
