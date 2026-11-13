export SInformation

"""
    SInformation <: HigherOrderInformationMeasure
    SInformation(; base = 2)

The S-information [Rosas2019](@cite): the total amount of information shared among
`n ≥ 2` variables, whether redundant or synergistic.

## Usage

- Use with [`association`](@ref) to compute the S-information of `n` input variables,
    `association(est, x₁, x₂, …, xₙ)`.

## Compatible estimators

- [`JointProbabilities`](@ref)
- [`EntropyDecomposition`](@ref)
- [`GaussianDecomposition`](@ref)

## Definition

```math
\\Sigma(X_1, \\ldots, X_n) = TC + DTC = \\sum_{i=1}^n I(X_i; X_{-i}),
```

where ``TC`` is the [`TotalCorrelation`](@ref), ``DTC`` is the
[`DualTotalCorrelation`](@ref), ``I`` is the Shannon mutual information and ``X_{-i}``
denotes all variables except ``X_i``.

See also: [`OInformation`](@ref).
"""
Base.@kwdef struct SInformation{B} <: HigherOrderInformationMeasure
    base::B = 2
end

function entropy_combination(::SInformation, H, V)
    return entropy_combination(TotalCorrelation(), H, V) +
        entropy_combination(DualTotalCorrelation(), H, V)
end
