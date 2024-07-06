using Test
using Random

@testset "PoczosSchneiderCMI" begin
    rng = Xoshiro(1234)
    x = rand(rng, 200)
    y = rand(rng, 200)
    z = rand(rng, 200)
    est = PoczosSchneiderCMI(CMIRenyiPoczos())
    @test association(est, x, y, z) isa Real
end