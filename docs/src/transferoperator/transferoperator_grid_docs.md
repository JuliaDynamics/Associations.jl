# Transfer operator approximation

## Grid estimator

For longer time series, we use a rectangular grid to discretize the embedding.
This approach is relatively fast, because no intersections between volumes
have to be explicitly computed.

### Documentation

```@docs
transferoperator_grid
```
