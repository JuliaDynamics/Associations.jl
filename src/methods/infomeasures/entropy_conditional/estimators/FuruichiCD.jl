# Todo: is the continuous Tsallis conditional entropy well defined and can be computed
# using a difference of Tsallis entropies? Then we could use EntropyEstimator.
export FuruichiCD

"""
    FuruichiCD <: ConditionalEntropyEstimator
    Furuichi(est::ProbabilitiesEstimator)

The `FuruichiCD` estimator (Furuichi, 2006)[^Furuichi2006] compute the discrete Tsallis
conditional entropy by first computing probabilities using `est`, then pluggin those
probabilities into the formula below.

## Description

The `FuruichiCD` estimator can use any [`ProbabilitiesEstimator`](@ref) to compute a plug-in
estimate of

```math
H_q^T(Y | X)
= - \\sum_{x, y} p(x, y)^q  \\log_q p(x|y)
= H_q^T(X, Y) - H_q^T(X)
```

where ``\\log_q`` is the `q`-logarithm (see Furuichi, 2006), and `q != 1`.

!!! warn "Common outcome space"
    X and Y must share outcome space for these estimates to make sense.

[^Furuichi2006]:
    Furuichi, S. (2006). Information theoretical properties of Tsallis
    entropies. Journal of Mathematical Physics, 47(2), 023302.
"""
struct FuruichiCD{P<:ProbabilitiesEstimator} <: ConditionalEntropyEstimator
    est::P
end

function entropy_conditional(e::Tsallis, est::FuruichiCD, x, y)
    entropy(e, est.est, Dataset(x, y)) - entropy(e, est.est, Dataset(x))
end
