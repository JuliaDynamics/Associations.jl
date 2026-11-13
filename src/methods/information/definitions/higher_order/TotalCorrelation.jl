export TotalCorrelation

"""
    TotalCorrelation <: HigherOrderInformationMeasure
    TotalCorrelation(; base = 2)

The total correlation [Watanabe1960](@cite), also called the multi-information:
the amount of information shared among `n ≥ 2` variables, counting each piece of
redundant information as many times as it is repeated.

## Usage

- Use with [`association`](@ref) to compute the total correlation of `n` input variables,
    `association(est, x₁, x₂, …, xₙ)`.

## Compatible estimators

- [`JointProbabilities`](@ref)
- [`EntropyDecomposition`](@ref)
- [`GaussianDecomposition`](@ref)

## Definition

```math
TC(X_1, \\ldots, X_n) = \\sum_{i=1}^n H(X_i) - H(X_1, \\ldots, X_n),
```

where ``H`` is the Shannon (joint) entropy. ``TC = 0`` if and only if the variables are
independent. For two variables, ``TC`` equals the mutual information
([`MIShannon`](@ref)).

See also: [`DualTotalCorrelation`](@ref), [`OInformation`](@ref), [`SInformation`](@ref).
"""
Base.@kwdef struct TotalCorrelation{B} <: HigherOrderInformationMeasure
    base::B = 2
end

function entropy_combination(::TotalCorrelation, H, V)
    return sum(H((i,)) for i in V) - H(V)
end
