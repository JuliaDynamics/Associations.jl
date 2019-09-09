# Cross map over multiple time series lengths

To perform a convergent cross map analysis as in [^1] one can apply the `crossmap` functions on time series of increasing length. The
`convergentcrossmap(driver, response, timeserieslengths; kwargs...)` function
does so by applying `crossmap` for each time series length in
`timeserieslengths`, where time windows always start at the first data point.


## Documentation

```@docs
CrossMappings.convergentcrossmap
```

[^1]:
    Sugihara, G., May, R., Ye, H., Hsieh, C. H., Deyle, E., Fogarty, M., & Munch, S. (2012). Detecting causality in complex ecosystems. Science. [https://doi.org/10.1126/science.1227079](https://doi.org/10.1126/science.1227079)