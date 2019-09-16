using DynamicalSystems

# Short time series of two coupled logistic maps (x drives y)
sys = logistic2_unidir(c_xy = 0.5)
npts = 40
x, y = columns(trajectory(sys, npts, Ttr = 200))

# average over a few different binning schemes, use a few different
# prediction lags and use the transfer operator grid estimator on the point
# cloud generated from the invariant measure over the triangulation.
binnings = [RectangularBinning(i) for i = 2:3]
estimator = TransferOperatorGrid()
ηs = 1:3

# Perform causality test, both from x to y and y to x
test = ApproximateSimplexIntersectionTest(binning = binnings, estimator = estimator, ηs = ηs)

@test causality(x, y, test) isa Vector{T} where {T}
