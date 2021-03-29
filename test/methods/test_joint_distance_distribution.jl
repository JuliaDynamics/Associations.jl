x, y = rand(1000), rand(1000)

import HypothesisTests: OneSampleTTest

@test joint_distance_distribution(x, y) isa Vector
@test joint_distance_distribution(OneSampleTTest, x, y) isa OneSampleTTest
