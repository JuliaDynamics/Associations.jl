import UncertainData: RandomSequences 

x = rand(100)
y = rand(100)

test_ccm = ConvergentCrossMappingTest(timeseries_lengths = [45, 50], n_reps = 10)
test_cm = CrossMappingTest(n_reps = 10)
test_vf = VisitationFrequencyTest(binning = RectangularBinning(5), ηs = 1:5)
test_tog = TransferOperatorGridTest(binning = RectangularBinning(5), ηs = 1:5)
test_jdd = JointDistanceDistributionTest()
test_jddt = JointDistanceDistributionTTest()
predtest = VisitationFrequencyTest(binning = RectangularBinning(5), ηs = -5:5)
test_pa = PredictiveAsymmetryTest(predictive_test = predtest)
test_npa = NormalisedPredictiveAsymmetryTest(predtest, f = 1.0)
test_NNmi = NearestNeighbourMITest(ηs = -5:5)

Nseq = 80
nchunks = 2
rtest_cm = RandomSequencesTest(test_cm, RandomSequences(nchunks, Nseq))
rtest_ccm = RandomSequencesTest(test_ccm, RandomSequences(nchunks, Nseq))
rtest_vf = RandomSequencesTest(test_vf, RandomSequences(nchunks, Nseq))
rtest_tog = RandomSequencesTest(test_tog, RandomSequences(nchunks, Nseq))
rtest_jdd = RandomSequencesTest(test_jdd, RandomSequences(nchunks, Nseq))
rtest_jddt = RandomSequencesTest(test_jddt, RandomSequences(nchunks, Nseq))
rtest_pa = RandomSequencesTest(test_pa, RandomSequences(nchunks, Nseq))
rtest_npa = RandomSequencesTest(test_npa, RandomSequences(nchunks, Nseq))
rtest_NNmi = RandomSequencesTest(test_NNmi, RandomSequences(nchunks, Nseq))

@test causality(test_cm, x, y) isa CausalAnalysis
@test causality(test_ccm, x, y) isa CausalAnalysis
@test causality(test_vf, x, y) isa CausalAnalysis
@test causality(test_tog, x, y) isa CausalAnalysis
@test causality(test_jdd, x, y) isa CausalAnalysis
@test causality(test_jddt, x, y) isa CausalAnalysis
@test causality(test_pa, x, y) isa CausalAnalysis
@test causality(test_npa, x, y) isa CausalAnalysis
@test causality(test_NNmi, x, y) isa CausalAnalysis