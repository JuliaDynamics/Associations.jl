export TsallisMIFuruichi

"""
    TsallisMIFuruichi <: MutualInformationEstimator
    TsallisMIFuruichi(est::ProbabilityEstimator)

The `TsallisMIFuruichi` estimator computes the discrete Tsallis mutual "information"
(called conditional mutual Tsallis *entropy* in Furuichi, 2006; see their
paper for reasoning on naming the method).

## Description

`TsallisMIFuruichi` is compatible with any [`ProbabilitiesEstimator`](@ref). You just need
to make sure that `est` is compatible with your input data.

Given a probabilities estimator, the `TsallisMIFuruichi` estimator first computes
probabilities using `est`, then plugs these probabilities into the formula below.

```math
I_q^T(X; Y) = H_q^T(X) - H_q^T(X | Y) = H_q^T(X) + H_q^T(Y) - H_q^T(X, Y)
```
"""
struct TsallisMIFuruichi{P <: ProbabilitiesEstimator} <: MutualInformationEstimator
    est::P
end

function mutualinfo(e::Tsallis, est::TsallisMIFuruichi, x, y)
    @assert e.q > 1
    entropy(e, est.est, Dataset(x)) +
        entropy(e, est.est, Dataset(y)) -
        entropy(e, est.est, Dataset(x, y))
end
