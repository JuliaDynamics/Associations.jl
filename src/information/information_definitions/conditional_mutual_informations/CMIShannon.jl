using ComplexityMeasures: Shannon
import ComplexityMeasures: log_with_base

export CMIShannon

"""
    CMIShannon <: ConditionalMutualInformation
    CMIShannon(; base = 2)

The Shannon conditional mutual information (CMI) ``I^S(X; Y | Z)``.

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise dependence.
- Use with [`condmutualinfo`](@ref) to compute the raw conditional mutual information.

## Supported definitions

Consider random variables ``X \\in \\mathbb{R}^{d_X}`` and
``Y \\in \\mathbb{R}^{d_Y}``, given ``Z \\in \\mathbb{R}^{d_Z}``. The Shannon
conditional mutual information is defined as

```math
\\begin{align*}
I(X; Y | Z)
&= H^S(X, Z) + H^S(Y, z) - H^S(X, Y, Z) - H^S(Z) \\\\
&= I^S(X; Y, Z) + I^S(X; Y)
\\end{align*},
```

where ``I^S(\\cdot; \\cdot)`` is the Shannon mutual information [`MIShannon`](@ref),
and ``H^S(\\cdot)`` is the [`Shannon`](@ref) entropy.

Differential Shannon CMI is obtained by replacing the entropies by
differential entropies.

See also: [`condmutualinfo`](@ref).
"""
Base.@kwdef struct CMIShannon{B} <: ConditionalMutualInformation
    base::B = 2
end

function information(definition::CMIShannon, pxyz::Probabilities{T, 3}) where T
    dx, dy, dz = size(pxyz)
    pxz = marginal(pxyz, dims = [1, 3])
    pyz = marginal(pxyz, dims = [2, 3])
    pz = marginal(pxyz, dims = 3)
    cmi = 0.0
    log0 = log_with_base(definition.base)
    for k in 1:dz
        pzₖ = pz[k]
        for j in 1:dy
            pyⱼzₖ = pyz[j, k]
            pyⱼzₖ > 0 || continue # leads to NaN
            for i in 1:dx
                pxᵢzₖ = pxz[i, k]
                pxᵢzₖ > 0 || continue # leads to NaN
                pxᵢyⱼzₖ = pxyz[i, j, k]
                inner = (pzₖ * pxᵢyⱼzₖ) / (pxᵢzₖ * pyⱼzₖ)
                if inner != 0.0
                    cmi += pxᵢyⱼzₖ * log0(inner)
                end
            end
        end
    end
    return cmi
end

# ------------------------------------------------
# Four-entropies decompostion of CMIShannon
# ------------------------------------------------
function information(est::EntropyDecomposition{<:CMIShannon, <:DifferentialInfoEstimator{<:Shannon}}, x, y, z)
    HXZ, HYZ, HXYZ, HZ = marginal_entropies_cmi4h_differential(est, x, y, z)
    cmi = HXZ + HYZ - HXYZ - HZ
    return cmi
end

function information(est::EntropyDecomposition{<:CMIShannon, <:DiscreteInfoEstimator{<:Shannon}}, x, y, z)
    HXZ, HYZ, HXYZ, HZ = marginal_entropies_cmi4h_discrete(est, x, y, z)
    cmi = HXZ + HYZ - HXYZ - HZ
    return cmi
end

# ---------------------------------------------------
# Two-mutual-information decomposition of CMIShannon 
# ---------------------------------------------------
function information(est::MIDecomposition{<:ConditionalMutualInformation, <:MutualInformationEstimator{<:MIShannon}}, x, y, z)
    MI_X_YZ, MI_X_Z = marginal_mutual_informations(est, x, y, z)
    cmi = MI_X_YZ - MI_X_Z
    return cmi
end

# We don't care if the estimated is mixed, discrete or handles both. The MI estimator 
# handles that internally.
function marginal_mutual_informations(est::MIDecomposition{<:ConditionalMutualInformation, <:MutualInformationEstimator{<:MIShannon}}, x, y, z)
    X = StateSpaceSet(x)
    Y = StateSpaceSet(y)
    Z = StateSpaceSet(z)
    YZ = StateSpaceSet(Y, Z)

    modified_est = estimator_with_overridden_parameters(est.definition, est.est)
    MI_X_YZ = information(modified_est, X, YZ)
    MI_X_Z = information(modified_est, X, Z)

    return MI_X_YZ, MI_X_Z
end

# ------------------------------------------------
# Pretty printing for decomposition estimators.
# ------------------------------------------------
function decomposition_string(definition::CMIShannon, 
        est::EntropyDecomposition{M, <:DiscreteInfoEstimator}) where {M}
    return "CMI_S(X, Y | Z) = H_S(X,Z) + H_S(Y,Z) - H_S(X,Y,Z) - H_S(Z)";
end

function decomposition_string(definition::CMIShannon, 
        est::EntropyDecomposition{M, <:DifferentialInfoEstimator}) where {M}
    return "CMI_S(X, Y | Z) = h_S(X,Z) + h_S(Y,Z) - h_S(X,Y,Z) - h_S(Z)";
end

function decomposition_string(definition::CMIShannon, 
    est::MIDecomposition{<:MIShannon, <:MutualInformationEstimator}) where {M}
    return "CMI_S(X, Y | Z) = I_S(X; Y, Z) + I_S(X; Z)"
end