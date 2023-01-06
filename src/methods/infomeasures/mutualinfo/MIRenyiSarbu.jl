export MIRenyi

"""
    MIRenyi <: MutualInformation
    MIRenyi(; base = 2, q = 1.5)

The discrete Rényi mutual information (see [`MRenyi`](@ref)) from Sarbu (2014)[^Sarbu2014].

## Description

Sarbu (2014) defines discrete Rényi mutual information as the
Rényi ``\\alpha``-divergence between the conditional joint probability mass function
``p(x, y)`` and the product of the conditional marginals, ``p(x) \\cdot p(y)``:

```math
I(X, Y; Z)^R_q =
\\dfrac{1}{q-1}
\\log \\left(
    \\sum{x \\in X}\\sum{y \\in Y}
    \\dfrac{p(x, y)^q}{\\left( p(x)\\cdot p(y) \\right)^{q-1}}
\\right)
```

!!! note "Only valid for fixed outcome spaces"
    This estimator is only well defined for [`ProbabilitiesEstimator`](@ref)s with fixed
    outcomes spaces. This is synonymous with having to estimate probabilities in the
    *joint* space, the marginalizing to obtain joint and marginal probability distributions
    with the same number of elements. This works out-of-the-box for some estimators
    like [`ValueHistogram`](@ref) with [`FixedRectangularBinning`](@ref), while for
    estimators such as [`SymbolicPermutation`](@ref), we provide specalized implementations.

[^Sarbu2014]: Sarbu, S. (2014, May). Rényi information transfer: Partial Rényi transfer
    entropy and partial Rényi mutual information. In 2014 IEEE International Conference
    on Acoustics, Speech and Signal Processing (ICASSP) (pp. 5666-5670). IEEE.

See also: [`mutualinfo`](@ref).
"""
struct MIRenyiSarbu{E <: Renyi} <: MutualInformation{E}
    e::E
    function MIRenyiSarbu(; q = 1.5, base = 2)
        e = Renyi(; q, base)
        new{typeof(e)}(e)
    end
end

function estimate(measure::MIRenyiSarbu, pxy::ContingencyMatrix{T, 2}) where {T}
    pxy = ContingencyMatrix(est, x, y)
    px = marginal_probs(pxy, 1)
    py = marginal_probs(pxy, 2)

    mi = 0.0
    for i in eachindex(px)
        for j in eachindex(py)
            pxyᵢⱼ = pxy[i, j]
            mi += pxyᵢⱼ^q / ((px[i] * py[j])^(q - 1))
        end
    end
    if mi == 0
        return 0.0
    else
        return (1 / (q - 1) * log(mi)) / log(e.base, ℯ)
    end
end
