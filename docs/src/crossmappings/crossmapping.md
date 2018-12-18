# Cross mappings

## Cross mapping algorithm

Given two data series,
the putative `driver` and the putative `response`, the the `crossmap(driver, response; kwargs...)` function computes how well a delay
embedding of `response` predicts scalar values of `driver`. This is [Sugihara et al. (2012)](http://science.sciencemag.org/content/early/2012/09/19/science.1227079)'s cross mapping algorithm.


To perform lagged CCM analysis ([Ye et al., 2015)](https://www.nature.com/articles/srep14750)), you can tune the `ν` parameter.

### Documentation

```@docs
CrossMappings.crossmap
```

## Convergent cross mapping (CCM)

To perform a convergent cross map analysis as in [Sugihara et al. (2012)](http://science.sciencemag.org/content/early/2012/09/19/science.1227079),
one can apply the `crossmap` functions on time series of increasing length. The
`convergentcrossmap(driver, response, timeserieslengths; kwargs...)` function
does so by applying `crossmap` for each time series length in
`timeserieslengths`, where time windows always start at the first data point.


### Documentation

```@docs
CrossMappings.convergentcrossmap
```


## Related software

- **[CauseMap.jl](https://github.com/cyrusmaher/CauseMap.jl)** also provides
    an implementation of the CCM algorithm, but this package has not been
    updated since 2015.
