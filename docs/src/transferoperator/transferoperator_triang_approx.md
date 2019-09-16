# Transfer operator estimators

## Approximate triangulation estimator

For short time series, the most reliable estimates of the transfer operator are
obtained by using a triangulation of the state space as our partition. This
approach is computationally costly because it has to compute N-dimensional
simplex intersections. This estimator computes approximate simplex intersections, which
is roughly an order of magnitude faster than computing exact intersections. This
estimator is suitable for datasets with a number of points in the order of hundreds.
For larger data sets, use a rectangular estimator.

The estimator returns a `TransferOperatorTriangulationApprox` instance.

### Documentation

```@docs
transferoperator_triangulation_approx(invariant_pts)
```
