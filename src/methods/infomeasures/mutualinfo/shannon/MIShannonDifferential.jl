export MIShannonDifferential
export ShannonH3Differential

"""
    ShannonH3Differential <: Definition

`ShannonH3Differential` is a directive used in combination with a multivariate
[`EntropyEstimator`](@ref) to compute the continuous Shannon mutual information between
``X \\in \\mathbb{R}^{d_X}`` and``Y \\in \\mathbb{R}^{d_Y}`` using the formula

```math
I^S(X; Y) = H^S(X) + H_q^S(Y) - H^S(X, Y),
```

where ``XY \\in \\mathbb{R}^{d_X + d_Y}`` is the joint space, and ``H^S(\\cdot)`` and
``H^S(\\cdot, \\cdot)`` are the marginal and joint Shannon entropies. Unless
otherwise specified, estimation is done in the most naive way possible:
compute probability mass functions separately for each space `X`, `Y` and `XY`, then
plug these probabilites into the respective entropy formulas.

## Supported estimators

- [`Kraskov`](@ref)
- [`KozachenkoLeonenko`](@ref)
- [`GaoNaive`](@ref) and [`GaoNaiveCorrected`](@ref)
- [`Goria`](@ref)
- [`Zhu`](@ref)
- [`ZhuSingh`](@ref)
- [`LeonenkoProzantoSavani`](@ref)
- [`Lord`](@ref)

See also: [`mutualinfo`](@ref).
"""
struct ShannonH3Differential <: Definition end

"""
    MIShannonDifferential <: MutualInformation
    MIShannonDifferential(e::Entropy = Shannon(), definition = ShannonH3Differential())
    MIShannonDifferential(; base::Real)

`MIShannonDifferential` is a directive to compute the differential Shannon mutual
information to base `e.b` between continuous variables ``X \\in \\mathbb{R}^{d_X}`` and
``Y \\in \\mathbb{R}^{d_Y}`` using the formula.

```math
I^S(X; Y) = h^S(X) + h_q^S(Y) - h^S(X, Y).
```

## Supported definitions

- [`ShannonH3Differential`](@ref).

See also: [`mutualinfo`](@ref).
"""
struct MIShannonDifferential{D <: Definition, E <: Renyi} <: MutualInformation
    e::E
    definition::D
    function MIShannonDifferential(; base = 2,
            definition::D = ShannonH3()) where {D}
            e = Shannon(; base)
        new{D, typeof(e)}(e, definition)
    end
end

function estimate(def::MIShannonDifferential{D, E}, est::EntropyEstimator, x, y) where {D, E}
    e = def.e
    X = Dataset(x)
    Y = Dataset(y)
    XY = Dataset(X, Y)
    return entropy(e, est, X) + entropy(e, est, Y) - entropy(e, est, XY)
end
