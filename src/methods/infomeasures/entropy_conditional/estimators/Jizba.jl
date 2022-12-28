export Jizba
# Todo: is the continuous Renyi conditional entropy well defined and can be computed
# using a difference of Renyi entropies? Then we could use DifferentialEntropyEstimator.
"""
    Jizba <: ConditionalDifferentialEntropyEstimator
    Jizba(est::ProbabilitiesEstimator)

Compute the discrete conditional Rényi entropy by first estimating probabilities using
`est`, then plugging them into the formula below.

## Description

The `Jizba` estimator can use any [`ProbabilitiesEstimator`](@ref) to compute a plug-in
estimate of

```math
R_q(Y | X)
= \\dfrac{1}{q - 1} \\log \\dfrac{\\sum_{x, y} P(x, y)^q}{\\sum_x P(X)^q}
= R_q(X, Y) - R_q(X).
```

[^Jizba2004]:
    P. Jizba and T. Arimitsu, “The world according to Rényi: Thermodynamics of
    multifractal systems,” Ann. Phys., vol. 312, pp. 17–59, 2004.
"""
struct Jizba{P <: ProbabilitiesEstimator} <: ConditionalDifferentialEntropyEstimator
    pest::P
end

function entropy_conditional(e::Renyi, est::Jizba, x, y)
    entropy(e, est.pest, Dataset(x, y)) - entropy(e, est.pest, Dataset(x))
end
