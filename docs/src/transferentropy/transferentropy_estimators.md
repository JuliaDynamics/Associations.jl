# Low-level transfer entropy estimators

For complete control over the estimation procedure, the analyst must create a [delay 
reconstruction](@ref custom_delay_reconstruction) from the input data, specify a [discretization scheme](@ref discretization) which can be either [rectangular](@ref partition_rectangular) or [triangulated](@ref partition_triangulated), and [map variables of the delay reconstruction to the correct 
marginals](@ref TEVars).

The package provides some convenience methods to compute TE directly from time series, see the [wrappers](@ref wrapper_TE) for
[regular TE](@ref wrapper_te_reg) and [conditional TE](@ref wrapper_te_cond). However, be absolutely sure that you understand what they do before applying them to real problems.

## [Estimators](@id te_estimators)

Valid estimator types for rectangular partitions are 

- `VisitationFrequency`. An implementation of the original TE estimator from Schreiber (2000)[^1]
- `TransferOperatorGrid`. An implementation of the transfer operator grid estimator from Diego et al. (2019)[^2]

## [General workflow for TE estimation](@id general_workflow_te)

The general workflow for estimating transfer entropy over rectangular partitions is as follows.

## [Rectangular partitions](@id te_estimator_rectangular)

To estimate transfer entropy over rectangular partitions, you would use the following method, providing
either a `VisitationFrequency` or `TransferOperatorGrid` instance to the `estimator` argument.

```@docs
transferentropy(pts, vars::TEVars, ϵ::RectangularBinning, estimator::TransferEntropyEstimator; b = 2)
```


## [Triangulated partitions](@id te_estimator_triangulated)

Estimators for computing TE on triangulated partitions, whose invariant distribution is obtained 
through the transfer operator, was also introduced in Diego et al. (2019)[^2].

```@docs
transferentropy(μ::AbstractTriangulationInvariantMeasure, vars::TEVars,
        binning_scheme::RectangularBinning;
        estimator = VisitationFrequency(), n::Int = 10000, b = 2)
```

[^1]:
    Schreiber, Thomas. "Measuring information transfer." Physical review letters 85.2 (2000): 461.
[^2]:
    Diego, David, Kristian Agasøster Haaga, and Bjarte Hannisdal. "Transfer entropy computation using the Perron-Frobenius operator." Physical Review E 99.4 (2019): 042212.
