# Cross mapping algorithms

## Cross map over multiple time series lengths

To perform a convergent cross map analysis as in [Sugihara et al. (2012)](http://science.sciencemag.org/content/early/2012/09/19/science.1227079),
one can apply the `crossmap` functions on time series of increasing length. The
`convergentcrossmap(driver, response, timeserieslengths; kwargs...)` function
does so by applying `crossmap` for each time series length in
`timeserieslengths`, where time windows always start at the first data point.


### Documentation

```@docs
CrossMappings.convergentcrossmap
```
