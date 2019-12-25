# [Transfer entropy estimation procedure](@id te_estimation_procedure)

## [Assigning marginals](@id te_assigning_marginals)

Assume we want to compute transfer entropy from a source variable ``S`` 
to a target variable ``T``, potentially conditioned on variable(s) ``C``.  
To control the transfer entropy, the analyst must 

1. Create a [generalised delay reconstruction](@ref custom_delay_reconstruction) from 
    the input data. In practice, this means that lagged instances of the time series 
    are organised in some data structure that can be sampled column-wise (e.g. a 
    [custom delay reconstruction](@ref custom_delay_reconstruction) or a `Dataset`).
2. Specify a [discretization scheme](@ref discretization), which can be either 
    [rectangular](@ref partition_rectangular) or [triangulated](@ref partition_triangulated).
3. Map variables of the generalised delay reconstruction to the correct marginals. In 
    other words, which columns of the input data should be grouped in which marginals? 
    There are four marginals that must be assigned:

    - The future of the target (``T_{f}``).
    - The present and past of the target variable (``T_{pp}``).
    - The present and past of the source variable (``S_{pp}``).
    - The present/past/future of any conditional variables (``C_{pp}``). If ``C_{pp}`` is non-empty,
        then conditional transfer entropy will be computed. Otherwise, non-conditioned transfer 
        entropy is computed.

The latter is achieved by creating a [`TEVars`](@ref) instance, which refers to the columns of the 
input data by their indices.

```@docs
TEVars(; )
```

## [Generalised delay reconstructions (convenience function)](@id te_embed_funcs)

```@docs
te_embed
```

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
