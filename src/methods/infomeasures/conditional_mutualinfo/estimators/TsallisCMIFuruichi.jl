export TsallisCMIFuruichi

"""
    TsallisCMIFuruichi <: MutualInformationEstimator
    TsallisCMIFuruichi(est::ProbabilityEstimator)

The `TsallisCMIFuruichi` estimator computes the discrete Tsallis conditional mutual
"information" (called conditional mutual Tsallis *entropy* in Furuichi, 2006; see their
paper for reasoning on naming the method).

## Description

`TsallisCMIFuruichi` is compatible with any [`ProbabilitiesEstimator`](@ref). You just need
to make sure that `est` is compatible with your input data.

Given a probabilities estimator, the `TsallisMIFuruichi` estimator first computes
probabilities using `est`, then plugs these probabilities into the formula below.

```math
I_q^T(X; Y | Z)
= H_q^T(X; Z) + H_q^T(X; Y, Z)
= I_q^T(X; Y, Z) - I_q^T(X; Z)
```

[^Furuichi2006]:
    Furuichi, S. (2006). Information theoretical properties of Tsallis
    entropies. Journal of Mathematical Physics, 47(2), 023302.
"""
struct TsallisCMIFuruichi{P <: ProbabilitiesEstimator} <: MutualInformationEstimator
    est::P
end

# TODO: use estimate interface with MI2
function cmi(e::Tsallis, est::TsallisCMIFuruichi, x, y, z)
    @assert e.q > 1
    X = Dataset(X)
    return mutualinfo(e, est.est, X, Dataset(y, z)) -
        mutualinfo(e, est.est, X, Dataset(z))
end
