using Test
using Random

@testset "JointProbabilities" begin
    rng = Xoshiro(1234)
    # Pre-discretized
    x = rand(rng, ["a", "b", "c"], 200)
    y = rand(rng, ["hello", "yoyo", "heyhey"], 200)
    z = rand(rng, ["a", "b"], 200)
    est = JointProbabilities(CMIShannon(), UniqueElements())
    @test association(est, x, y, z) >= 0.0
    @test association(est, x, y, z) isa Real
end