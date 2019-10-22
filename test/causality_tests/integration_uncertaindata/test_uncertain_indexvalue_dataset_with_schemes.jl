using UncertainData

N = 50
x_uncertain = [UncertainValue(Normal, x, rand(Uniform(0.1, 0.8))) for x in rand(N)]
y_uncertain = [UncertainValue(Normal, y, rand(Uniform(0.1, 0.8))) for y in rand(N)]
x = UncertainValueDataset(x_uncertain)
y = UncertainValueDataset(y_uncertain)

time_uncertain = [UncertainValue(Normal, i, 1) for i = 1:length(x)];
time_certain = [CertainValue(i) for i = 1:length(x)];
timeinds_x = UncertainIndexDataset(time_uncertain)
timeinds_y = UncertainIndexDataset(time_certain)

X = UncertainIndexValueDataset(timeinds_x, x)
Y = UncertainIndexValueDataset(timeinds_y, y);

η_max = 3
ηs = -η_max:η_max
k, l, m = 1, 1, 1
n_bins = ceil(Int, N^(1/(k + l + m + 1)))

test_vf = VisitationFrequencyTest(k = 1, l = 1, m = 1, τ = 2, binning = RectangularBinning(n_bins), ηs = 1:η_max)
test_tog = TransferOperatorGridTest(k = 1, l = 1, m = 1, τ = 2, binning = RectangularBinning(n_bins), ηs = 1:η_max)
test_te = TransferOperatorGridTest(k = 1, l = 1, m = 1, τ = 2, binning = RectangularBinning(n_bins), ηs = ηs)
test_pa = PredictiveAsymmetryTest(predictive_test = test_te)

test_ccm = ConvergentCrossMappingTest(timeseries_lengths = [45, 50], n_reps = 20)
test_cm = CrossMappingTest(n_reps = 10)
test_jdd = JointDistanceDistributionTest()
test_jddt = JointDistanceDistributionTTest()

# Truncate each of the indices for x at 0.8 their standard deviation around the mean
constraints_x_inds = TruncateStd(0.8)

# Truncate each of the indices for y at 1.5 their standard deviation around the mean
constraints_y_inds = TruncateStd(1.5)

# Truncate each of the values of x at the 20th percentile range
constraints_x_vals = [TruncateQuantiles(0.4, 0.6) for i = 1:N];

# Truncate each of the values of x at the 80th percentile range
constraints_y_vals = [TruncateQuantiles(0.1, 0.9) for i = 1:N];
cs_x = (constraints_x_inds, constraints_x_vals)
cs_y = (constraints_y_inds, constraints_y_vals)
grid = RegularGrid(1, N, 1)

seq_type = StrictlyIncreasing()
n_draws = 10
resampling_constraints = ConstrainedIndexValueResampling(n_draws, cs_x, cs_y);
resampling_sequential = SequentialResampling(seq_type)
resampling_sequential_intp = SequentialInterpolatedResampling(seq_type, grid)

@test causality(X, Y, test_pa, resampling_constraints) isa Vector{Vector{T}} where T <: Real
@test causality(X, Y, test_pa, resampling_constraints, resampling_sequential) isa Vector{Vector{T}} where T <: Real
@test causality(X, Y, test_pa, resampling_constraints, resampling_sequential_intp) isa Vector{Vector{T}} where T <: Real

@test causality(X, Y, test_vf, resampling_constraints) isa Vector{Vector{T}} where T <: Real
@test causality(X, Y, test_vf, resampling_constraints, resampling_sequential) isa Vector{Vector{T}} where T <: Real
@test causality(X, Y, test_vf, resampling_constraints, resampling_sequential_intp) isa Vector{Vector{T}} where T <: Real

@test causality(X, Y, test_tog, resampling_constraints) isa Vector{Vector{T}} where T <: Real
@test causality(X, Y, test_tog, resampling_constraints, resampling_sequential) isa Vector{Vector{T}} where T <: Real
@test causality(X, Y, test_tog, resampling_constraints, resampling_sequential_intp) isa Vector{Vector{T}} where T <: Real

@test causality(X, Y, test_cm, resampling_constraints) isa Vector{Vector{T}} where T <: Real
@test causality(X, Y, test_cm, resampling_constraints, resampling_sequential) isa Vector{Vector{T}} where T <: Real
@test causality(X, Y, test_cm, resampling_constraints, resampling_sequential_intp) isa Vector{Vector{T}} where T <: Real

@test causality(X, Y, test_ccm, resampling_constraints) isa Vector{Vector{Vector{T}}} where T <: Real
@test causality(X, Y, test_ccm, resampling_constraints, resampling_sequential) isa Vector{Vector{Vector{T}}} where T <: Real
@test causality(X, Y, test_ccm, resampling_constraints, resampling_sequential_intp) isa Vector{Vector{Vector{T}}} where T <: Real

@test causality(X, Y, test_jdd, resampling_constraints) isa Vector{Vector{T}} where T <: Real
@test causality(X, Y, test_jdd, resampling_constraints, resampling_sequential) isa Vector{Vector{T}} where T <: Real
@test causality(X, Y, test_jdd, resampling_constraints, resampling_sequential_intp) isa Vector{Vector{T}} where T <: Real

@test causality(X, Y, test_jddt, resampling_constraints) isa Vector{OneSampleTTest}
@test causality(X, Y, test_jddt, resampling_constraints, resampling_sequential) isa Vector{OneSampleTTest}
@test causality(X, Y, test_jddt, resampling_constraints, resampling_sequential_intp) isa Vector{OneSampleTTest}