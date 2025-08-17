using ComplexityMeasures: Tsallis

export MITsallisMartin

"""
    MITsallisMartin <: BivariateInformationMeasure
    MITsallisMartin(; base = 2, q = 1.5)

The discrete Tsallis mutual information from [Martin2004](@citet).

## Usage

- Use with [`association`](@ref) to compute the raw Tsallis-Martin mutual information from input data
    using of of the estimators listed below.
- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise dependence using
    the Tsallis-Martin mutual information.

## Compatible estimators

- [`JointProbabilities`](@ref)
- [`EntropyDecomposition`](@ref)

## Description

Martin et al.'s Tsallis mutual information between variables ``X \\in \\mathbb{R}^{d_X}``
and ``Y \\in \\mathbb{R}^{d_Y}`` is defined as

```math
I_{\\text{Martin}}^T(X, Y, q) := H_q^T(X) + H_q^T(Y) - (1 - q) H_q^T(X) H_q^T(Y) - H_q(X, Y),
```

where ``H^S(\\cdot)`` and ``H^S(\\cdot, \\cdot)`` are the marginal and joint Shannon
entropies, and `q` is the [`Tsallis`](@extref ComplexityMeasures)-parameter.

## Estimation

- [Example 1](@ref example_MITsallisMartin_JointProbabilities_UniqueElements): [`JointProbabilities`](@ref) with [`UniqueElements`](@extref ComplexityMeasures) outcome space.
- [Example 2](@ref example_MITsallisMartin_EntropyDecomposition_LeonenkoProzantoSavani): [`EntropyDecomposition`](@ref) with [`LeonenkoProzantoSavani`](@extref ComplexityMeasures) estimator.
- [Example 3](@ref example_MITsallisMartin_EntropyDecomposition_OrdinalPatterns): [`EntropyDecomposition`](@ref) with [`OrdinalPatterns`](@extref ComplexityMeasures) outcome space.
"""
Base.@kwdef struct MITsallisMartin{B,Q} <: MutualInformation
    base::B = 2
    q::Q = 1.5
end


# ----------------------------------------------------------------
# Estimation methods
# ----------------------------------------------------------------
function association(est::JointProbabilities{<:MITsallisMartin}, x, y)
    probs = probabilities(est.discretization, x, y)
    return association(est.definition, probs)
end

# This is definition 3 in Martin et al. (2004), but with pᵢ replaced by the joint
# distribution and qᵢ replaced by the product of the marginal distributions.
function association(definition::MITsallisMartin, pxy::Probabilities{T,2}) where T
    (; base, q) = definition
    # TODO: return MIShannon if q = 1? otherwise, we don't need `base`.
    q != 1 || throw(ArgumentError("`MITsallisMartin` for q=$(q) not defined."))
    px = marginal(pxy, dims=1)
    py = marginal(pxy, dims=2)

    mi = 0.0
    for (i, pxᵢ) in enumerate(px.p)
        for (j, pyⱼ) in enumerate(py.p)
            pxyᵢⱼ = pxy[i, j]
            mi += pxyᵢⱼ^q / (pxᵢ^(q - 1) * pyⱼ^(q - 1))
        end
    end
    f = 1 / (q - 1)
    return f * (1 - mi)
end

function association(est::EntropyDecomposition{<:MITsallisMartin,<:DifferentialInfoEstimator{<:Tsallis}}, x, y)
    HX, HY, HXY = marginal_entropies_mi3h_differential(est, x, y)
    q = est.definition.q
    mi = HX + HY - (1 - q) * HX * HY - HXY
    return mi
end

function association(est::EntropyDecomposition{<:MITsallisMartin,<:DiscreteInfoEstimator{<:Tsallis}}, x, y)
    HX, HY, HXY = marginal_entropies_mi3h_discrete(est, x, y)
    q = est.definition.q
    mi = HX + HY - (1 - q) * HX * HY - HXY
    return mi
end

# ------------------------------------------------
# Pretty printing for decomposition estimators.
# ------------------------------------------------
function decomposition_string(definition::MITsallisMartin, est::DiscreteInfoEstimator{<:Tsallis})
    return "MI_S(X, Y) = H_T(X) + H_T(Y) - (1 - q)*H_T(X)*H_T(Y) - H_T(X, Y)"
end

function decomposition_string(definition::MITsallisMartin, est::DifferentialInfoEstimator{<:Tsallis})
    return "MI_S(X, Y) = h_T(X) + h_T(Y) - (1 - q)*h_T(X)*H_T(Y) - h_T(X, Y)"
end