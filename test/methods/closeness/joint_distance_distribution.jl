using Random
import HypothesisTests: OneSampleTTest

rng = MersenneTwister(123456)
x, y = rand(rng, 1000), rand(rng, 1000)

@test jdd(x, y) isa Vector
@test jdd(OneSampleTTest, x, y) isa OneSampleTTest


# v2.X and upwards
@test association(JointDistanceDistribution(), x, y) isa Vector

@test_logs (:warn, "Convenience function `jdd` is deprecated. Use `association(JointDistanceDistribution(; kwargs...), x, y)` instead.") jdd(x, y)
@test_logs (:warn, "Convenience function `jdd` is deprecated. Use `association(JointDistanceDistribution(; kwargs...), x, y)` instead.") jdd(JointDistanceDistribution(), x, y)
@test_logs (:warn, "jdd(::OneSampleTTest, x, y; kwargs...) is deprecated. Instead, do `measure = JointDistanceDistribution(); independence(JointDistanceDistributionTest(measure), x, y)`.") jdd(OneSampleTTest, x, y)