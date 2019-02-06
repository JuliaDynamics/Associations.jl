# Cross mapping algorithms

## Cross mapping

Given two data series,
the putative `driver` and the putative `response`, the the `crossmap(driver, response; kwargs...)` function computes how well a delay
embedding of `response` predicts scalar values of `driver`. This is [Sugihara et al. (2012)](http://science.sciencemag.org/content/early/2012/09/19/science.1227079)'s cross mapping algorithm.


To perform lagged CCM analysis ([Ye et al., 2015)](https://www.nature.com/articles/srep14750)), you can tune the `ν` parameter.

### Documentation

```@docs
CrossMappings.crossmap
```
