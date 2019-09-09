import HypothesisTests: OneSampleTTest

x, y = rand(1000), rand(1000)
@test causality(x, y, JointDistanceDistributionTest()) isa Vector
@test causality(x, y, JointDistanceDistributionTTest()) isa OneSampleTTest