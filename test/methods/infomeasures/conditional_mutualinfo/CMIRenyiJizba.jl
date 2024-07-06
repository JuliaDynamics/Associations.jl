@testset "JointProbabilities" begin
    rng = Xoshiro(1234)
    # Pre-discretized
    x = rand(rng, ["a", "b", "c"], 200)
    y = rand(rng, ["hello", "yoyo", "heyhey"], 200)
    z = rand(rng, ["a", "b"], 200)
    est = JointProbabilities(CMIRenyiJizba(), UniqueElements())
    @test association(est, x, y, z) >= 0.0
    @test association(est, x, y, z) isa Real
end

@testset "EntropyDecomposition" begin
    rng = Xoshiro(1234)
    # Pre-discretized
    x = rand(rng, 200)
    y = rand(rng, 200)
    z = rand(rng, 200)
    @testset "Discrete" begin
        est_disc = EntropyDecomposition(CMIRenyiJizba(), PlugIn(Renyi(q = 2)), OrdinalPatterns(m=3))
        @test association(est_disc, x, y, z) >= 0.0
        @test association(est_disc, x, y, z) isa Real
    end

    @testset "Differential" begin
        est_diff = EntropyDecomposition(CMIRenyiJizba(), LeonenkoProzantoSavani(Renyi(q = 2)))
        @test association(est_disc, x, y, z) >= 0.0
        @test association(est_disc, x, y, z) isa Real
    end
end