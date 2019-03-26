using HypothesisTests

x = rand(500)
y = rand(500)

jdd_test = joint_distance_distribution(x, y)
@test jdd_test isa OneSampleTTest