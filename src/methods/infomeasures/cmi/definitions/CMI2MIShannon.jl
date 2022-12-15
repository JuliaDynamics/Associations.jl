export CMI2MI

"""
    CMI2MIShannon <: ConditionalMutualInformationDefinition
    CMI2MIShannon()

A directive to compute conditional mutual information (CMI) as a difference of mutual
information terms.

```math
I(X; Y | Z) = I^S(X; Y, Z) + I^S(X; Y)
```

where ``H^S(\\cdot)`` is the Shannon entropy of the argument. The formula applies both
for discrete and continuous CMI.
```
"""
struct CMI2MI <: ConditionalMutualInformationDefinition end

function estimate(def::CMI2MI, measure::CMIShannon, est::MutualInformationEstimator, x, y, z)
    X = Dataset(x)
    Y = Dataset(y)
    Z = Dataset(z)
    YZ = Dataset(Y, Z)
    return mutualinfo(MIShannon(), est, X, YZ) +
        mutualinfo(MIShannon(), est, X, Y)
end

estimate(measure::CMIShannon,
    est::MutualInformationEstimator, x, y, z) = estimate(CMI2MI(), measure, est, x, y, z)
