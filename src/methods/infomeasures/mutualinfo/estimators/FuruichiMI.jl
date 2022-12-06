export FuruichiMI

"""
    FuruichiMI <: MutualInformationEstimator
    FuruichiMI(est::ProbabilityEstimator)

The `FuruichiMI` estimator computes the discrete Tsallis mutual "information" (called
mutual Tsallis *entropy* in Furuichi, 2006). It does so by first computing probabilities
using `est`, then plugging these probabilities into the formula below.

## Description

```math
I_q^T(X; Y) = H_q^T(X) - H_q^T(X | Y) = H_q^T(X) + H_q^T(Y) - H_q^T(X, Y)
```
"""
struct FuruichiMI{P <: ProbabilitiesEstimator} <: MutualInformationEstimator
    est::P
end

function mutualinfo(e::Tsallis, est::FuruichiMI, x, y)
    @assert e.q > 1
    entropy(e, est.est, Dataset(x)) +
        entropy(e, est.est, Dataset(y)) -
        entropy(e, est.est, Dataset(x, y))
end
