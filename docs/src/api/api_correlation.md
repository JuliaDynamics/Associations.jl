# [Correlation API](@id correlation_api)

## Pearson correlation

```@docs
PearsonCorrelation
pearson_correlation
```

## Partial correlation

```@docs
PartialCorrelation
partial_correlation
```

## Distance correlation

```@docs
DistanceCorrelation
distance_correlation
```

## Examples


### Distance correlation

Let's do a comparison with an example from the
[*energy*](https://cran.r-project.org/web/packages/energy/index.html) R-package.

```@example
using CausalityTools
x = -1.0:0.1:1.0 |> collect
y = map(xᵢ -> xᵢ^3 - 2xᵢ^2 - 3, x)
z = map(yᵢ -> yᵢ^2 - 2yᵢ, y)
dcov = distance_correlation(x, y)
round(dcov, digits = 3) == 0.673 # should be true
```