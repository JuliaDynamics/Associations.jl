using Test
using CausalityTools
using Random
rng = MersenneTwister(1234)

# ----------------------------------------------------------------
# Joint probabilities estimation
# ----------------------------------------------------------------
# all MI measures can be computed from the joint pmf
definitions = [MIShannon(), MIRenyiJizba(), MIRenyiSarbu(), MITsallisFuruichi(), MITsallisMartin()]

@testset "JointProbabilities with $(typeof(def).name.name)" for def in definitions
    x, y = rand(rng, 100), rand(rng, 100)
    X, Y = StateSpaceSet(rand(rng, 100, 2)), StateSpaceSet(rand(rng, 100, 2))
    
    est = JointProbabilities(def, ValueBinning(2), RelativeAmount())
    # The estimation of probabilities is decoupled from the estimation of the mutual info.
    # We could in principle use any probabilities estimator here, but we default to `RelativeAmount`.
    @test association(est, x, Y) isa Real;
    @test association(est, X, y) isa Real;
    @test association(est, x, y) isa Real;
    @test association(est, X, Y) isa Real;
end

# ----------------------------------------------------------------
# Decomposition estimation
# ----------------------------------------------------------------

# The following measures can be estimated using an entropy decomposition
defs = [MIShannon(), MITsallisMartin(), MITsallisFuruichi()]
ests = [PlugIn(Shannon()), PlugIn(Tsallis(q = 1.5)), PlugIn(Tsallis(q = 1.5))]
@testset "EntropyDecomposition with $(typeof(def).name.name)" for (def, est) in zip(defs, ests)
    x, y = rand(rng, 100), rand(rng, 100)
    X, Y = StateSpaceSet(rand(rng, 100, 2)), StateSpaceSet(rand(rng, 100, 2))
    est = EntropyDecomposition(def, est, OrdinalPatterns(m=2), RelativeAmount())
    #@show est
    # The estimation of probabilities is decoupled from the estimation of the mutual info.
    # We could in principle use any probabilities estimator here, but we default to `RelativeAmount`.
    @test association(est, x, Y) isa Real
    @test association(est, X, y) isa Real
    @test association(est, x, y) isa Real
    @test association(est, X, Y) isa Real
end

# The following measures cannot be decomposed into entropies and should throw errors
definitions = [MIRenyiSarbu()]
@testset "EntropyDecomposition with $(typeof(def).name.name)" for def in definitions
    x, y = rand(rng, 100), rand(rng, 100)
    X, Y = StateSpaceSet(rand(rng, 100, 2)), StateSpaceSet(rand(rng, 100, 2))
    
    est_diff = EntropyDecomposition(def, Kraskov(k=3))
    est_disc = EntropyDecomposition(def, PlugIn(Shannon()), OrdinalPatterns(m=2), RelativeAmount())

    @test_throws ArgumentError association(est_diff, x, Y)
    @test_throws ArgumentError association(est_diff, X, y)
    @test_throws ArgumentError association(est_diff, x, y)
    @test_throws ArgumentError association(est_diff, X, Y)

    @test_throws ArgumentError association(est_disc, x, Y)
    @test_throws ArgumentError association(est_disc, X, y)
    @test_throws ArgumentError association(est_disc, x, y)
    @test_throws ArgumentError association(est_disc, X, Y)
end
