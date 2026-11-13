export DualTotalCorrelation

"""
    DualTotalCorrelation <: HigherOrderInformationMeasure
    DualTotalCorrelation(; base = 2)

The dual total correlation [Han1978](@cite), also called the binding information:
the part of the joint entropy of `n ≥ 2` variables that is shared by at least two of
them, counting each piece of shared information once.

## Usage

- Use with [`association`](@ref) to compute the dual total correlation of `n` input
    variables, `association(est, x₁, x₂, …, xₙ)`.

## Compatible estimators

- [`JointProbabilities`](@ref)
- [`EntropyDecomposition`](@ref)
- [`GaussianDecomposition`](@ref)

## Definition

```math
DTC(X_1, \\ldots, X_n) = H(X_1, \\ldots, X_n) - \\sum_{i=1}^n H(X_i \\mid X_{-i})
= \\sum_{i=1}^n H(X_{-i}) - (n - 1) H(X_1, \\ldots, X_n),
```

where ``H`` is the Shannon (joint) entropy and ``X_{-i}`` denotes all variables except
``X_i``. ``DTC = 0`` if and only if the variables are independent. For two variables,
``DTC`` equals the mutual information ([`MIShannon`](@ref)).

See also: [`TotalCorrelation`](@ref), [`OInformation`](@ref), [`SInformation`](@ref).
"""
Base.@kwdef struct DualTotalCorrelation{B} <: HigherOrderInformationMeasure
    base::B = 2
end

function entropy_combination(::DualTotalCorrelation, H, V)
    n = length(V)
    return sum(H(setdiff(V, i)) for i in V) - (n - 1) * H(V)
end
