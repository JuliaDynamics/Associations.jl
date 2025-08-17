using ComplexityMeasures: Renyi

export MIRenyiSarbu

"""
    MIRenyiSarbu <: BivariateInformationMeasure
    MIRenyiSarbu(; base = 2, q = 1.5)

The discrete Rényi mutual information from [Sarbu2014](@citet).

## Usage

- Use with [`association`](@ref) to compute the raw Rényi-Sarbu mutual information from input data
    using of of the estimators listed below.
- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise dependence using
    the Rényi-Sarbu mutual information.

## Compatible estimators

- [`JointProbabilities`](@ref).

## Description

Sarbu (2014) defines discrete Rényi mutual information as the
Rényi ``\\alpha``-divergence between the conditional joint probability mass function
``p(x, y)`` and the product of the conditional marginals, ``p(x) \\cdot p(y)``:

```math
I(X, Y)^R_q =
\\dfrac{1}{q-1}
\\log \\left(
    \\sum_{x \\in X, y \\in Y}
    \\dfrac{p(x, y)^q}{\\left( p(x)\\cdot p(y) \\right)^{q-1}}
\\right)
```

## Estimation

- [Example 1](@ref example_MIRenyiSarbu_JointProbabilities_UniqueElements): [`JointProbabilities`](@ref) with [`UniqueElements`](@extref ComplexityMeasures) for categorical data.
- [Example 2](@ref example_MIRenyiSarbu_JointProbabilities_CosineSimilarityBinning): [`JointProbabilities`](@ref) with [`CosineSimilarityBinning`](@ref) for numerical data.
"""
Base.@kwdef struct MIRenyiSarbu{B,Q} <: MutualInformation
    base::B = 2
    q::Q = 1.5
end

# ----------------------------------------------------------------
# Estimation methods
# ----------------------------------------------------------------
function association(est::JointProbabilities{<:MIRenyiSarbu}, x, y)
    probs = probabilities(est.discretization, x, y)
    return association(est.definition, probs)
end

function association(definition::MIRenyiSarbu, pxy::Probabilities{T,2}) where T
    (; base, q) = definition

    px = marginal(pxy, dims=1)
    py = marginal(pxy, dims=2)

    mi = 0.0
    for i in eachindex(px.p)
        for j in eachindex(py.p)
            pxyᵢⱼ = pxy[i, j]
            mi += pxyᵢⱼ^q / ((px[i] * py[j])^(q - 1))
        end
    end
    if mi == 0
        return 0.0
    else
        return _convert_logunit(1 / (q - 1) * log(mi), ℯ, base)
    end
end

function association(est::EntropyDecomposition{<:MIRenyiSarbu}, x, y)
    throw(ArgumentError("MIRenyiSarbu not implemented for $(typeof(est).name.name)"))
end
