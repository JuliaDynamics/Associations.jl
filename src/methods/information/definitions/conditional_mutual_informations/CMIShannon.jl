using ComplexityMeasures: Shannon
import ComplexityMeasures: log_with_base

export CMIShannon

"""
    CMIShannon <: ConditionalMutualInformation
    CMIShannon(; base = 2)

The Shannon conditional mutual information (CMI) ``I^S(X; Y | Z)``.

## Usage

- Use with [`association`](@ref) to compute the raw Shannon conditional mutual information
    using of of the estimators listed below.
- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise conditional 
    independence using the Shannon conditional mutual information.

## Compatible estimators

- [`JointProbabilities`](@ref)
- [`EntropyDecomposition`](@ref)
- [`MIDecomposition`](@ref)
- [`FPVP`](@ref)
- [`MesnerShalizi`](@ref)
- [`Rahimzamani`](@ref)
- [`PoczosSchneiderCMI`](@ref)
- [`GaussianCMI`](@ref)

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

## Estimation

- [Example 1](@ref example_CMIShannon_EntropyDecomposition_Kraskov): 
    [`EntropyDecomposition`](@ref) with [`Kraskov`](@extref ComplexityMeasures) estimator.
- [Example 2](@ref example_CMIShannon_EntropyDecomposition_ValueBinning):
    [`EntropyDecomposition`](@ref) with [`ValueBinning`](@extref ComplexityMeasures) estimator.
- [Example 3](@ref example_CMIShannon_MIDecomposition_KSG1): 
    [`MIDecomposition`](@ref) with [`KraskovStögbauerGrassberger2`](@extref ComplexityMeasures) estimator.
"""
Base.@kwdef struct CMIShannon{B} <: ConditionalMutualInformation
    base::B = 2
end

# ----------------------------------------------------------------
# Estimation methods
# ----------------------------------------------------------------
function association(est::JointProbabilities{<:CMIShannon}, x, y, z)
    probs = probabilities(est.discretization, x, y, z)
    return association(est.definition, probs)
end

function association(definition::CMIShannon, pxyz::Probabilities{T,3}) where T
    dx, dy, dz = size(pxyz)
    pxz = marginal(pxyz, dims=[1, 3])
    pyz = marginal(pxyz, dims=[2, 3])
    pz = marginal(pxyz, dims=3)
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
function association(est::EntropyDecomposition{<:CMIShannon,<:DifferentialInfoEstimator{<:Shannon}}, x, y, z)
    HXZ, HYZ, HXYZ, HZ = marginal_entropies_cmi4h_differential(est, x, y, z)
    cmi = HXZ + HYZ - HXYZ - HZ
    return cmi
end

function association(est::EntropyDecomposition{<:CMIShannon,<:DiscreteInfoEstimator{<:Shannon}}, x, y, z)
    HXZ, HYZ, HXYZ, HZ = marginal_entropies_cmi4h_discrete(est, x, y, z)
    cmi = HXZ + HYZ - HXYZ - HZ
    return cmi
end

# ---------------------------------------------------
# Two-mutual-information decomposition of CMIShannon 
# ---------------------------------------------------
function association(est::MIDecomposition{<:ConditionalMutualInformation,<:MutualInformationEstimator{<:MIShannon}}, x, y, z)
    MI_X_YZ, MI_X_Z = marginal_mutual_informations(est, x, y, z)
    cmi = MI_X_YZ - MI_X_Z
    return cmi
end

# We don't care if the estimated is mixed, discrete or handles both. The MI estimator 
# handles that internally.
function marginal_mutual_informations(est::MIDecomposition{<:ConditionalMutualInformation,<:MutualInformationEstimator{<:MIShannon}}, x, y, z)
    X = StateSpaceSet(x)
    Y = StateSpaceSet(y)
    Z = StateSpaceSet(z)
    YZ = StateSpaceSet(Y, Z)

    modified_est = estimator_with_overridden_parameters(est.definition, est.est)
    MI_X_YZ = association(modified_est, X, YZ)
    MI_X_Z = association(modified_est, X, Z)

    return MI_X_YZ, MI_X_Z
end

# ---------------------------------
# Avoid some common errors
# ---------------------------------
function verify_decomposition_entropy_type(definition::CMIShannon, est::INFO_ESTS)
    if !(est.definition isa Shannon)
        T = typeof(est.definition).name.name
        msg = "Can't decompose CMIShannon into a combination of $T entropies. Please provide a `Shannon` entropy estimator instead."
        throw(ArgumentError(msg))
    end
end


# ------------------------------------------------
# Pretty printing for decomposition estimators.
# ------------------------------------------------
function decomposition_string(
    definition::CMIShannon,
    est::EntropyDecomposition{<:CMIShannon,<:DiscreteInfoEstimator{<:Shannon}}
)
    return "Iₛ(X, Y | Z) = Hₛ(X,Z) + Hₛ(Y,Z) - Hₛ(X,Y,Z) - Hₛ(Z)"
end

function decomposition_string(
    definition::CMIShannon,
    est::EntropyDecomposition{<:CMIShannon,<:DifferentialInfoEstimator{<:Shannon}}
)
    return "Iₛ(X, Y | Z) = hₛ(X,Z) + hₛ(Y,Z) - hₛ(X,Y,Z) - hₛ(Z)"
end

function decomposition_string(
    definition::CMIShannon,
    est::MIDecomposition{<:CMIShannon,<:MutualInformationEstimator}
)
    return "Iₛ(X, Y | Z) = Iₛ(X; Y, Z) + Iₛ(X; Z)"
end