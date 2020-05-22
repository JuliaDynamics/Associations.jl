# Cross mapping

Given two data series,
the putative `driver` and the putative `response`, the the `crossmap(driver, response; kwargs...)` function computes how well a delay
embedding of `response` predicts scalar values of `driver`. This is the original cross mapping algorithm from Sugihara et al. [^1].


To perform lagged CCM analysis [^2] as Ye et al., you can tune the `ν` parameter.

## Documentation

```@docs
CrossMappings.crossmap
```

[^1]:
    Sugihara, G., May, R., Ye, H., Hsieh, C. H., Deyle, E., Fogarty, M., & Munch, S. (2012). Detecting causality in complex ecosystems. Science. [https://doi.org/10.1126/science.1227079](https://doi.org/10.1126/science.1227079)
[^2]:
    Ye, H., Deyle, E. R., Gilarranz, L. J., & Sugihara, G. (2015). Distinguishing time-delayed causal interactions using convergent cross mapping. Scientific Reports. [https://doi.org/10.1038/srep14750](https://doi.org/10.1038/srep14750)
