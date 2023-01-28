export CETsallisAbe

"""
    CETsallisAbe <: ConditionalEntropy
    CETsallisAbe(; base = 2, q = 1.5)

Abe & Rajagopal (2001)'s discrete Tsallis conditional entropy measure.

## Definition

Abe & Rajagopal's Tsallis conditional entropy between discrete random variables
``X`` and ``Y`` with finite ranges ``\\mathcal{X}`` and ``\\mathcal{Y}`` is defined as

```math
H_q^{T_A}(X | Y) = \\dfrac{H_q^T(X, Y) - H_q^T(Y)}{1 + (1-q)H_q^T(Y)},
```

where ``H_q^T(\\cdot)`` and ``H_q^T(\\cdot, \\cdot)`` is the [`Tsallis`](@ref)
entropy and the joint Tsallis entropy.

[^Abe2001]:
    Abe, S., & Rajagopal, A. K. (2001). Nonadditive conditional entropy and its
    significance for local realism. Physica A: Statistical Mechanics and its Applications,
    289(1-2), 157-164.
"""
struct CETsallisAbe{E} <: ConditionalEntropy
    e::E
    function CETsallisAbe(; q = 1.5, base = 2)
        e = MLEntropy(Tsallis(; q, base))
        new{typeof(e)}(e)
    end
end

function estimate(measure::CETsallisAbe, pxy::ContingencyMatrix{T, 2}) where {T}
    e = measure.e.definition
    Nx, Ny = size(pxy)
    base, q = e.base, e.q

    py = probabilities(pxy, 2)
    # Definition 7 in Abe & Rajagopal (2001)
    hjoint = 1 / (1 - q) * (sum(pxy .^ 2) - 1)

    # The marginal Tsallis entropy for the second variable
    hy = entropy(Tsallis(; q, base), py)

    # Equation 13 in Abe & Rajagopal (2001)
    ce = (hjoint - hy) / (1 + (1 - q)*hy)

    if q == 1 # if shannon, normalize
        return ce / log(base, ℯ)
    else
        return ce
    end
end

function estimate(measure::CETsallisAbe, est::ProbabilitiesEstimator, x, y)
    e = measure.e.definition
    q, base = e.q, e.base

    HY, HXY = marginal_entropies_ce2h(measure, est, x, y)
    ce = (HXY - HY) / (1 + (1 - q)*HY)
    if q == 1 # if shannon, normalize
        return ce / log(base, ℯ)
    else
        return ce
    end
end

function estimate(measure::CETsallisAbe, est::DifferentialEntropyEstimator, x, y)
    throw(ArgumentError("CETsallisAbe not implemented for $(typeof(est))"))
end
