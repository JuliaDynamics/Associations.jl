export FuruichiCMI

"""
    FuruichiCMI <: MutualInformationEstimator
    FuruichiCMI(est::ProbabilityEstimator)

The `FuruichiCMI` estimator computes the discrete Tsallis conditional mutual "information"
(called mutual Tsallis *entropy* in Furuichi, 2006). It does so by first computing probabilities
using `est`, then plugging these probabilities into the formula below.

## Description

```math
I_q^T(X; Y | Z)
= H_q^T(X; Z) + H_q^T(X; Y, Z)
= I_q^T(X; Y, Z) - I_q^T(X; Z)
```

"""
struct FuruichiCMI{P <: ProbabilitiesEstimator} <: MutualInformationEstimator
    est::P
end

# TODO: use estimate interface with MI2
function cmi(e::Tsallis, est::FuruichiCMI, x, y, z)
    @assert e.q > 1
    X = Dataset(X)
    return mutualinfo(e, est.est, X, Dataset(y, z)) -
        mutualinfo(e, est.est, X, Dataset(z))
end
