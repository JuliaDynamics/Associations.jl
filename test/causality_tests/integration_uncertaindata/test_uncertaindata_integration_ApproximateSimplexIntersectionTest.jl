using DynamicalSystems

# Short time series of two coupled logistic maps (x drives y)
sys = logistic2_unidir(c_xy = 0.5)
npts = 40
x, y = columns(trajectory(sys, npts, Ttr = 200))

# Add some uncertainties to the time series
uvals_x = [UncertainValue(Normal, X, abs(rand(Normal(0, 0.05)))) for X in x]
uvals_y = [UncertainValue(Normal, Y, abs(rand(Normal(0, 0.05)))) for Y in y]
UVX = UncertainValueDataset(uvals_x)
UVY = UncertainValueDataset(uvals_y)

# average over a few different binning schemes, use a few different
# prediction lags and use the transfer operator grid estimator on the point
# cloud generated from the invariant measure over the triangulation.
binnings = [RectangularBinning(i) for i = 3:4]
estimator = TransferOperatorGrid()
ηs = [1]

# Perform causality test, both from x to y and y to x
test = ApproximateSimplexIntersectionTest(binning = binnings, estimator = estimator, ηs = ηs)
@test causality(uvals_x, y, test) isa Vector{T} where {T}
@test causality(uvals_x, uvals_y, test) isa Vector{T} where {T}
@test causality(x, uvals_y, test) isa Vector{T} where {T}
@test causality(UVX, y, test) isa Vector{T} where {T}
@test causality(UVX, UVY, test) isa Vector{T} where {T}
@test causality(x, UVY, test) isa Vector{T} where {T}
