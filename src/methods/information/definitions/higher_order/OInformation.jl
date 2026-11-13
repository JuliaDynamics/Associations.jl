export OInformation

"""
    OInformation <: HigherOrderInformationMeasure
    OInformation(; base = 2)

The O-information [Rosas2019](@cite): a signed measure of whether the information
shared among `n` variables is dominated by redundancy (``Ω > 0``) or by synergy
(``Ω < 0``).

## Usage

- Use with [`association`](@ref) to compute the O-information of `n` input variables,
    `association(est, x₁, x₂, …, xₙ)`.
- Use with [`oinformation_gradients`](@ref) to find how much each variable (or pair of
    variables) contributes to the O-information.

## Compatible estimators

- [`JointProbabilities`](@ref)
- [`EntropyDecomposition`](@ref)
- [`GaussianDecomposition`](@ref)

## Definition

```math
\\Omega(X_1, \\ldots, X_n) = TC - DTC
= (n - 2) H(X_1, \\ldots, X_n) + \\sum_{i=1}^n \\left[ H(X_i) - H(X_{-i}) \\right],
```

where ``TC`` is the [`TotalCorrelation`](@ref), ``DTC`` is the
[`DualTotalCorrelation`](@ref), ``H`` is the Shannon (joint) entropy and ``X_{-i}``
denotes all variables except ``X_i``.

For example, ``n`` copies of one fair bit give ``Ω = n - 2`` bits (pure redundancy), and
``n - 1`` independent fair bits plus their parity give ``Ω = 2 - n`` bits (pure synergy).
For ``n = 2``, ``Ω = 0``. For ``n = 3``, ``Ω`` equals the co-information
``I(X_1; X_2) - I(X_1; X_2 \\mid X_3)``.

See also: [`SInformation`](@ref).
"""
Base.@kwdef struct OInformation{B} <: HigherOrderInformationMeasure
    base::B = 2
end

function entropy_combination(::OInformation, H, V)
    n = length(V)
    return (n - 2) * H(V) + sum(H((i,)) - H(setdiff(V, i)) for i in V)
end
