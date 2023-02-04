using HypothesisTests: OneSampleTTest
x, y = rand(100), rand(100)
test =  JointDistanceDistributionTest()
@test test isa JointDistanceDistributionTest
@test independence(test, x, y) isa OneSampleTTest
