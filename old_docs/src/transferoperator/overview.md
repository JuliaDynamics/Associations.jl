# Transfer operator estimation

There are currently three different transfer operator estimators implemented:

- [`TransferOperatorEstimatorTriangulationExact`](transferoperator_triang_exact.md) approximates
    the transfer operator from a triangulation using exact simplex intersections.
- [`transferoperator_triangulation_approx`](transferoperator_triang_approx.md)
    approximates the transfer operator from a triangulation using approximate simplex intersections.
- [`TransferOperatorEstimatorRectangularBinning`](transferoperator_rectangular_binning.md)
    approximates the transfer operator from a rectangular binning.

These estimators return instances of the following types.

- [`TransferOperatorTriangulationExact`](composite_types/triang_exact.md)
- [`TransferOperatorTriangulationApprox`](composite_types/triang_approx.md)
- [`TransferOperatorRectangularBinning`](composite_types/rectangular_binning.md)
