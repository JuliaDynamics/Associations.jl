# Invariant measure estimation

When the transfer operator has been precomputed for a state space discretization, we can
easily derive an invariant measure (probability density) over the states of the system by
calling `invariantmeasure` on the transfer operator.

This works both for transfer operators estimated from triangulations[^1] and from
rectangular partitions [^2].

## Documentation

```@docs
invariantmeasure
```

[^1]:
    Using the [`transferoperator_triangulation_exact`](../transferoperator/transferoperator_triang_exact.md) and [`transferoperator_triangulation_approx`](../transferoperator/transferoperator_triang_approx.md) estimators.

[^2]:
    Using the [`TransferOperatorEstimatorRectangularBinning`](../transferoperator/transferoperator_rectangular_binning.md) estimator.
