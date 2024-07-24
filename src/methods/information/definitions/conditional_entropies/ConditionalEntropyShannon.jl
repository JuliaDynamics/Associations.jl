using ComplexityMeasures: Shannon
import ComplexityMeasures: log_with_base

export ConditionalEntropyShannon

"""
    ConditionalEntropyShannon <: ConditionalEntropy
    ConditionalEntropyShannon(; base = 2)

The [`Shannon`](@ref) conditional entropy measure.

## Usage 

- Use with [`association`](@ref) to compute the Shannon conditional entropy between 
    two variables.

## Compatible estimators

- [`JointProbabilities`](@ref)

## Discrete definition

### Sum formulation

The conditional entropy between discrete random variables
``X`` and ``Y`` with finite ranges ``\\mathcal{X}`` and ``\\mathcal{Y}`` is defined as

```math
H^{S}(X | Y) = -\\sum_{x \\in \\mathcal{X}, y \\in \\mathcal{Y}} p(x, y) \\log(p(x | y)).
```

This is the definition used when calling [`association`](@ref) with a
[`JointProbababilities`](@ref) estimator.

### Two-entropies formulation

Equivalently, the following differenConditionalEntropy of entropies hold

```math
H^S(X | Y) = H^S(X, Y) - H^S(Y),
```

where ``H^S(\\cdot)`` and ``H^S(\\cdot | \\cdot)`` are the [`Shannon`](@ref) entropy and
Shannon joint entropy, respectively. This is the definition used when calling
[`association`](@ref) with a [`ProbabilitiesEstimator`](@ref).

## Differential definition

The differential conditional Shannon entropy is analogously defined as

```math
H^S(X | Y) = h^S(X, Y) - h^S(Y),
```

where ``h^S(\\cdot)`` and ``h^S(\\cdot | \\cdot)`` are the [`Shannon`](@ref)
differential entropy and Shannon joint differential entropy, respectively. This is the
definition used when calling [`association`](@ref) with a
[`DifferentialEntropyEstimator`](@extref).

## Estimation

- [Example 1](@ref example_ConditionalEntropyShannon_analytical): Analytical example from Cover & Thomas's book.
- [Example 2](@ref example_ConditionalEntropyShannon_JointProbabilities_CodifyVariables_UniqueElements): 
    [`JointProbabilities`](@ref) estimator with[`CodifyVariables`](@ref) discretization and 
    [`UniqueElements`](@ref) outcome space on categorical data.
- [Example 3](@ref example_ConditionalEntropyShannon_JointProbabilities_CodifyPoints_UniqueElementsEncoding): 
    [`JointProbabilities`](@ref) estimator with [`CodifyPoints`](@ref) discretization and [`UniqueElementsEncoding`](@ref)
    encoding of points on numerical data.
"""
Base.@kwdef struct ConditionalEntropyShannon{B} <: ConditionalEntropy
    base::B = 2
end

# ----------------------------------------------------------------
# Estimation methods
# ----------------------------------------------------------------
function association(est::JointProbabilities{<:ConditionalEntropyShannon}, inputs...)
    probs = probabilities(est.discretization, inputs...)
    return association(est.definition, probs)
end

function association(definition::ConditionalEntropyShannon, pxy::Probabilities{T, 2}) where {T}
    base = definition.base
    Nx, Ny = size(pxy)
    py = marginal(pxy, dims = 2)

    ce = 0.0
    log0 = log_with_base(base)
    for j in 1:Ny
        pyⱼ = py[j]
        for i in 1:Nx
            pxyᵢⱼ = pxy[i, j]
            if pxyᵢⱼ != 0.0
                ce += pxyᵢⱼ * log0(pxyᵢⱼ / pyⱼ)
            end
        end
    end
    return -ce
end

