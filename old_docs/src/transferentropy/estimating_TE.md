# [Transfer entropy estimation procedure](@id te_estimation_procedure)

## [Computing transfer entropy from rectangular partitions](@id te_estimator_rectangular)

```@docs
transferentropy(::Any, ::TEVars, ::RectangularBinning, ::BinningTransferEntropyEstimator)
```

## [Computing transfer entropy by counting nearest neighbors](@id te_estimator_nn)

```@docs
transferentropy(::Any, ::TEVars, ::NearestNeighbourMI)
```

## [Computing transfer entropy from triangulated partitions](@id te_estimator_triang)

Estimators for computing transfer entropy on triangulated partitions, whose invariant distribution is obtained through the transfer operator, was introduced in Diego et al. (2019)[^2].

```@docs
transferentropy(μ::AbstractTriangulationInvariantMeasure, vars::TEVars, binning_scheme::RectangularBinning)
```

[^2]:
    Diego, David, Kristian Agasøster Haaga, and Bjarte Hannisdal. "Transfer entropy computation using the Perron-Frobenius operator." Physical Review E 99.4 (2019): 042212.
