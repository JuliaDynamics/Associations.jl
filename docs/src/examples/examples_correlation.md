# Correlation examples

## Distance correlation

Let's do a comparison with an example from the *energy* R-package.

```@docs
x = -1.0:0.1:1.0 |> collect
y = map(xᵢ -> xᵢ^3 - 2xᵢ^2 - 3, x)
z = map(yᵢ -> yᵢ^2 - 2yᵢ, y)
dcov = distance_correlation(x, y)
@test round(dcov, digits = 3) == 0.673
```