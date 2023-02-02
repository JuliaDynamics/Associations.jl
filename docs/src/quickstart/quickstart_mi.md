# [Mutual information](@id quickstart_mutualinfo)

## [`MIShannon`](@ref) (differential)

The differential Shannon mutual information ([`MIShannon`](@ref)) can be estimated
using a dedicated mutual information estimator like [`KraskovStögbauerGrassberger2`](@ref).
These estimators typically apply some form of bias correction.

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
mutualinfo(KSG2(k = 5), x, y)
```

We can also estimate [`MIShannon`](@ref) by naively applying a [`DifferentialEntropyEstimator`](@ref), which doesn't apply any bias correction.

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
mutualinfo(Kraskov(k = 3), x, y)
```

## [`MIShannon`](@ref) (discrete, numerical)

A [`ValueHistogram`](@ref) estimator can be used to bin the data and compute
discrete Shannon mutual information.

```@example mi_demonstration
using CausalityTools

# Use the H3-estimation method with a discrete visitation frequency based 
# probabilities estimator over a fixed grid covering the range of the data,
# which is on [0, 1].
est = ValueHistogram(FixedRectangularBinning(0, 1, 5))
mutualinfo(est, x, y)
```

If you need access to the estimated joint probability mass function,
use a [`ContingencyMatrix`](@ref). This is slower, but convenient if you need to
investigate the probabilities manually.

```@example mi_demonstration
using CausalityTools
c = contingency_matrix(est, x, y)
mutualinfo(c)
```

## [`MIShannon`](@ref) (discrete, categorical)

The [`ContingencyMatrix`](@ref) approach can also be used with categorical data.
For example, let's compare the Shannon mutual information between the preferences
of a population sample with regards to different foods.

```@example mi_demonstration
using CausalityTools
n = 1000
preferences = rand(["neutral", "like it", "hate it"], n);
random_foods = rand(["water", "flour", "bananas", "booze", "potatoes", "beans", "soup"], n)
biased_foods = map(preferences) do preference
    if cmp(preference, "neutral") == 1
        return rand(["water", "flour"])
    elseif cmp(preference, "like it") == 1
        return rand(["bananas", "booze"])
    else
        return rand(["potatoes", "beans", "soup"])
    end
end

c_biased = contingency_matrix(preferences, biased_foods) 
c_random = contingency_matrix(preferences, random_foods) 
mutualinfo(c_biased), mutualinfo(c_random)
```
