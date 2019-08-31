x, y = rand(1000), rand(1000)

@test joint_distance_distribution(x, y) isa OneSampleTTest
