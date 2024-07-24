
# ----------------------------------------------------------------
# Joint probabilities estimation
# ----------------------------------------------------------------
definitions = [CMIShannon(), CMIRenyiSarbu(), CMITsallisPapapetrou()]

@testset "JointProbabilities with $(typeof(def).name.name)" for def in definitions
    x, y, z = rand(rng, 100), rand(rng, 100), rand(rng, 100)
    X, Y, Z = StateSpaceSet(rand(rng, 100, 2)), 
        StateSpaceSet(rand(rng, 100, 2)), 
        StateSpaceSet(rand(rng, 100, 2))
    
    est = JointProbabilities(def, ValueBinning(2), RelativeAmount())
    # The estimation of probabilities is decoupled from the estimation of the mutual info.
    # We could in principle use any probabilities estimator here, but we default to `RelativeAmount`.
    @test association(est, x, Y, z) isa Real;
    @test association(est, X, y, z) isa Real;
    @test association(est, x, y, z) isa Real;
    @test association(est, X, Y, Z) isa Real;
end

# Not defined for joint probabilities estimator
defs = [CMIRenyiJizba()]