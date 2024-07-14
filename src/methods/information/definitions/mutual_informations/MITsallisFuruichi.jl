using ComplexityMeasures: Tsallis

export MITsallisFuruichi
"""
    MITsallisFuruichi <: BivariateInformationMeasure
    MITsallisFuruichi(; base = 2, q = 1.5)

The discrete Tsallis mutual information from Furuichi (2006)[Furuichi2006](@cite), which
in that paper is called the *mutual entropy*.

## Usage

- Use with [`association`](@ref) to compute the raw Tsallis-Furuichi mutual information from input data
    using of of the estimators listed below.
- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise dependence using
    the Tsallis-Furuichi mutual information.

## Compatible estimators

- [`JointProbabilities`](@ref)
- [`EntropyDecomposition`](@ref)

## Description

Furuichi's Tsallis mutual entropy between variables ``X \\in \\mathbb{R}^{d_X}`` and
``Y \\in \\mathbb{R}^{d_Y}`` is defined as

```math
I_q^T(X; Y) = H_q^T(X) - H_q^T(X | Y) = H_q^T(X) + H_q^T(Y) - H_q^T(X, Y),
```

where ``H^T(\\cdot)`` and ``H^T(\\cdot, \\cdot)`` are the marginal and joint Tsallis
entropies, and `q` is the [`Tsallis`](@ref)-parameter.
```
"""
Base.@kwdef struct MITsallisFuruichi{B, Q} <: MutualInformation
    base::B = 2
    q::Q = 1.5
end

# ----------------------------------------------------------------
# Estimation methods
# ----------------------------------------------------------------
function association(est::JointProbabilities{<:MITsallisFuruichi}, x, y)
    probs = probabilities(est.discretization, x, y)
    return association(est.definition, probs)
end

function association(definition::MITsallisFuruichi, pxy::Probabilities{T, 2}) where T
    (; base, q) = definition

    px = marginal(pxy, dims = 1)
    py = marginal(pxy, dims = 2)

    mi = 0.0
    for i in eachindex(px.p)
        for j in eachindex(py.p)
            pxyᵢⱼ = pxy[i, j]
            mi += pxyᵢⱼ^q / (px[i]^(q - 1) * py[j]^(q - 1))
        end
    end
    mi = (1 / (q - 1) * (1 - mi) / (1-q))
    return _convert_logunit(mi, ℯ, base)
end



function association(est::EntropyDecomposition{<:MITsallisFuruichi, <:DifferentialInfoEstimator{<:Tsallis}}, x, y)
    HX, HY, HXY = marginal_entropies_mi3h(est, x, y)
    mi = HX + HY - HXY
    return mi
end

function association(est::EntropyDecomposition{<:MITsallisFuruichi, <:DiscreteInfoEstimator{<:Tsallis}}, x, y)
    HX, HY, HXY = marginal_entropies_mi3h_discrete(est, x, y)
    mi = HX + HY - HXY
    return mi
end


# ------------------------------------------------
# Pretty printing for decomposition estimators.
# ------------------------------------------------
function decomposition_string(definition::MITsallisFuruichi, est::DiscreteInfoEstimator{<:Tsallis})
    return "MI_TF(X, Y) = H_T(X) + H_T(Y) - H_T(X, Y)";
end

function decomposition_string(definition::MITsallisFuruichi, est::DifferentialInfoEstimator{<:Tsallis})
    return "MI_TF(X, Y) = h_T(X) + h_T(Y) - h_T(X, Y)";
end