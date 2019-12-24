using DynamicalSystems 
using UncertainData
using Distributions
using Test 

############################################################
# Define a system and generate some uncertain time series
############################################################
system = ar1_unidir(c_xy = 0.7)
n_points =  300
orbit = trajectory(system, n_points, Ttr = 1000)

# Add uncertainties to the time series values
x_uncertain = [UncertainValue(Normal, x, rand(Uniform(0.1, 0.8))) for x in orbit[:, 1]]
y_uncertain = [UncertainValue(Normal, y, rand(Uniform(0.1, 0.8))) for y in orbit[:, 2]]
x = UncertainValueDataset(x_uncertain)
y = UncertainValueDataset(y_uncertain)

# Add uncertainties to only the time indices for. Take the 
# time indices for y as certain.
time_uncertain = [UncertainValue(Normal, i, 1) for i = 1:length(x)];
time_certain = [CertainValue(i) for i = 1:length(x)];
timeinds_x = UncertainIndexDataset(time_uncertain)
timeinds_y = UncertainIndexDataset(time_certain)

X = UncertainIndexValueDataset(timeinds_x, x)
Y = UncertainIndexValueDataset(timeinds_y, y)

############
# Run tests
############
η_max = 6
ηs = -η_max:η_max
k, l, m = 1, 1, 1
n_bins = ceil(Int, n_points^(1/(k + l + m + 1)))

test_vf = VisitationFrequencyTest(k = 1, l = 1, m = 1, τ = 2, binning = RectangularBinning(n_bins), ηs = 1:η_max)
test_tog = TransferOperatorGridTest(k = 1, l = 1, m = 1, τ = 2, binning = RectangularBinning(n_bins), ηs = 1:η_max)
test_te = TransferOperatorGridTest(k = 1, l = 1, m = 1, τ = 2, binning = RectangularBinning(n_bins), ηs = ηs)
test_pa = PredictiveAsymmetryTest(predictive_test = test_te)

test_ccm = ConvergentCrossMappingTest(timeseries_lengths = [45, 50], n_reps = 20)
test_cm = CrossMappingTest(n_reps = 10)
test_jdd = JointDistanceDistributionTest()
test_jddt = JointDistanceDistributionTTest()
# Sample strictly increasing realisations of the indices, interpolate,
# then perform the causality test on the interpolated data.
tstep = 2
grid = RegularGrid(0, length(X), tstep)

@test causality(X, Y, test_pa, StrictlyIncreasing(), grid) isa Vector{T} where T <: Real
@test causality(X, Y, test_ccm, StrictlyIncreasing(), grid) isa Vector{Vector{T}} where T <: Real
@test causality(X, Y, test_cm, StrictlyIncreasing(), grid) isa Vector{T} where T <: Real
@test causality(X, Y, test_vf, StrictlyIncreasing(), grid) isa Vector{T} where T <: Real
@test causality(X, Y, test_tog, StrictlyIncreasing(), grid) isa Vector{T} where T <: Real
@test causality(X, Y, test_jdd, StrictlyIncreasing(), grid) isa Vector{T} where T <: Real
@test causality(X, Y, test_jddt, StrictlyIncreasing(), grid) isa OneSampleTTest

