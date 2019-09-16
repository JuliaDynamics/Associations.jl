import HypothesisTests: OneSampleTTest

jtest = JointDistanceDistributionTTest()

@test causality(x, y, jtest) isa OneSampleTTest