import HypothesisTests: OneSampleTTest

# Example system
sys = logistic2_unidir(c_xy = 1.0)

# Run test on time series with 200 observations, treating the first 
# variable as the source and the second variable as the target
setup = DiscreteSystemSetup(source = 1, target = 2, n_pts = 200)

# Causality tests
test_ccm = ConvergentCrossMappingTest(timeseries_lengths = [45, 50], n_reps = 20)
test_cm = CrossMappingTest(n_reps = 10)
test_vf = VisitationFrequencyTest(binning = RectangularBinning(5), ηs = 1:5)
test_tog = TransferOperatorGridTest(binning = RectangularBinning(5), ηs = 1:5)
test_jdd = JointDistanceDistributionTest()
test_jddt = JointDistanceDistributionTTest()
predtest = VisitationFrequencyTest(binning = RectangularBinning(5), ηs = -5:5)
test_pa = PredictiveAsymmetryTest(predictive_test = predtest)

# Run test 
@test causality(sys, setup, test_ccm) isa Vector{Vector{T}} where {T <: Real}
@test causality(sys, setup, test_cm) isa Vector{T} where {T <: Real}
@test causality(sys, setup, test_vf) isa Vector{T} where {T <: Real}
@test causality(sys, setup, test_tog) isa Vector{T} where {T <: Real}
@test causality(sys, setup, test_jdd) isa Vector{T} where {T <: Real}
@test causality(sys, setup, test_jddt) isa OneSampleTTest
