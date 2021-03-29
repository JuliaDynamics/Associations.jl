x, y = rand(1000), rand(1000)

import HypothesisTests: OneSampleTTest

@test jdd(x, y) isa Vector
@test jdd(OneSampleTTest, x, y) isa OneSampleTTest
