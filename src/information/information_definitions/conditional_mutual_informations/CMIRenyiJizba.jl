using ComplexityMeasures: Renyi

export CMIRenyiJizba

"""
    CMIRenyiJizba <: ConditionalMutualInformation

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
"""
struct CMIRenyiJizba{E <: Renyi} <: ConditionalMutualInformation
    e::E
    function CMIRenyiJizba(; base = 2, q = 1.5)
        e = Renyi(; base, q)
        new{typeof(e)}(e)
    end
    function CMIRenyiJizba(e::E) where E <: Renyi
        new{E}(e)
    end
end

function information(est::JointProbabilities{<:CMIRenyiJizba}, x, y, z)
    throw(ArgumentError("CMIRenyiJizba not implemented for $(typeof(est).name.name)"))
end

# --------------------------------------------------------------
# Conditional mutual information through entropy decomposition
# --------------------------------------------------------------
function information(est::DifferentialDecomposition{<:CMIRenyiJizba, <:DifferentialInfoEstimator{<:Renyi}}, x, y, z)
    HXZ, HYZ, HXYZ, HZ = marginal_entropies_cmi4h_differential(est, x, y, z)
    cmi = HXZ + HYZ - HXYZ - HZ
    return cmi
end

function information(est::DiscreteDecomposition{<:CMIRenyiJizba, <:DiscreteInfoEstimator{<:Renyi}}, x, y, z)
    HXZ, HYZ, HXYZ, HZ = marginal_entropies_cmi4h_discrete(est, x, y, z)
    cmi = HXZ + HYZ - HXYZ - HZ
    return cmi
end

# ------------------------------------------------
# Pretty printing for decomposition estimators.
# ------------------------------------------------
function decomposition_string(definition::CMIRenyiJizba, est::DiscreteInfoEstimator)
    return "CMI_RJ(X, Y) = H_R(X,Z) + H_R(Y,Z) - H_R(X,Y,Z) - H_R(Z)";
end

function decomposition_string(definition::CMIRenyiJizba, est::DifferentialInfoEstimator)
    return "CMI_RJ(X, Y) = h_R(X,Z) + h_R(Y,Z) - h_R(X,Y,Z) - h_R(Z)";
end