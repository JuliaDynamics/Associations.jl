# Transfer operator estimators

## Exact triangulation estimator

For short time series, the most reliable estimates of the transfer operator are
obtained by using a triangulation of the state space as our partition. This
approach is computationally costly because it has to compute N-dimensional
simplex intersections. This estimator computes exact simplex intersections, which
is very time consuming. It should only be used for small data sets.

The estimator returns a `TransferOperatorTriangulationExact` instance.

### Documentation

```@docs
transferoperator_triangulation_exact(invariant_pts)
```
