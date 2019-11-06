using UncertainData 
import Interpolations: Linear, Flat, OnGrid
using Interpolations

# Define the InterpolateAndBin instance
intp_bin = InterpolateAndBin(mean, 0:10:100, Linear(), 0:0.01:100, Flat(OnGrid()))

ib_test = InterpolateBinTest(CrossMappingTest(), intp_bin)
@test ib_test isa InterpolateBinTest
@test ib_test.n == 1

ib_test = InterpolateBinTest(CrossMappingTest(), intp_bin, 10)
@test ib_test isa InterpolateBinTest
@test ib_test.n == 10

# Some test data
N = 100
sys = ar1_unidir(c_xy = 0.2)
X, Y = example_uncertain_indexvalue_datasets(sys, N, (1, 2),
    d_xval = Uniform(0.001, 0.2), d_yval = Uniform(0.001, 0.2));

# Define a causality test 
k, l, m = 1, 1, 1 # embedding parameters
n_subdivisions = floor(Int, N^(1/(k + l + m + 1)))
state_space_binning = RectangularBinning(n_subdivisions)

te_test = VisitationFrequencyTest(k = k, l = l, m = m,
            binning = state_space_binning, 
            ηs = ηs, b = 2) # use base-2 logarithms
pa_test = PredictiveAsymmetryTest(predictive_test = te_test)

# Define interpolation grid over the range of available index values
tmin = max(minimum(mean.(X.indices)), minimum(mean.(Y.indices)))
tmax = max(maximum(mean.(X.indices)), maximum(mean.(Y.indices)))
intp_grid = tmin:0.01:tmax

# Define binning grid
left_bin_edges = tmin:1.5:tmax

# Define the InterpolateAndBin instance
intp_bin = InterpolateAndBin(mean, left_bin_edges, Linear(), intp_grid, Flat(OnGrid()))

# Define interpolate-and-bin test
ib_test = InterpolateBinTest(pa_test, intp_bin, 2)

res_xy = causality(X, Y, ib_test)

@test res_xy isa Vector{Vector{T}} where T
@test length(res_xy) == ib_test.n
