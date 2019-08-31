x, y = rand(1000), rand(1000)

import HypothesisTests: OneSampleTTest

@test joint_distance_distribution(x, y) isa OneSampleTTest
