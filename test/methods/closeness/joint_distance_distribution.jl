using Random
import HypothesisTests: OneSampleTTest

rng = MersenneTwister(123456)
x, y = rand(rng, 1000), rand(rng, 1000)

@test jdd(x, y) isa Vector
@test jdd(OneSampleTTest, x, y) isa OneSampleTTest


# v2.X and upwards
@test association(JointDistanceDistribution(), x, y) isa Vector
