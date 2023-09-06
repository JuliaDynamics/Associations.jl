# TODO

- When using `ValueBinning`, provide a `FixedRectangularBinning` with the same
    dimensions as the joint StateSpaceSet. Or, provide a `RectangularBinning`, and first encode the joint data, then take the marginals of that.
- When using `TransferOperator`, what to do?

# Generic estimators
- `RankTransformed`: a directive used with any estimator to first rank-transform the input data.