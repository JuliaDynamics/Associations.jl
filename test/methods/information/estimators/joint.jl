using Associations
using Test
using Random
rng = Xoshiro(1234)

# ------------------------------------------------------------
# API test 
# ------------------------------------------------------------
x = randn(rng, 100)
y = x .+ randn(rng, 100)
z = y .+ randn(rng, 100)

os = [
    OrdinalPatterns(m=2),
    Dispersion(),
    ValueBinning(2),
]

bivariate_symmetric_measures = [
    MIShannon(),
    MITsallisMartin(),
    MITsallisFuruichi(),
    MIRenyiSarbu(),
    MIRenyiJizba(),
    JointEntropyRenyi(),
    JointEntropyTsallis(),
    JointEntropyShannon(),
]

@testset "JointProbabilities estimator with $(typeof(m).name.name)" for m in bivariate_symmetric_measures
    @testset "CodifyVariables with $(typeof(o).name.name)" for o in os
        est = JointProbabilities(m, CodifyVariables(o))
        est_xy = association(est, x, y)
        est_yx = association(est, y, x)
        
        @test est_xy isa Real
        @test est_yx isa Real
        @test est_xy ≈ est_yx # symmetry
    end
end;

bivariate_asymmetric_measures = [
    ConditionalEntropyShannon(),
    ConditionalEntropyTsallisAbe(),
    ConditionalEntropyTsallisFuruichi(),
    HellingerDistance(),
    KLDivergence(),
    RenyiDivergence(),
    VariationDistance(),
]

@testset "JointProbabilities estimator with $(typeof(m).name.name)" for m in bivariate_asymmetric_measures
    @testset "CodifyVariables with $(typeof(o).name.name)" for o in os
        est = JointProbabilities(m, CodifyVariables(o))
        @test association(est, x, y) isa Real
    end
end;

trivariate_asymmetric_measures = [
    CMIShannon(),
    CMITsallisPapapetrou(),
    CMIRenyiSarbu(),
]

@testset "JointProbabilities estimator with $(typeof(m).name.name)" for m in trivariate_asymmetric_measures
    @testset "CodifyVariables with $(typeof(o).name.name)" for o in os
        est = JointProbabilities(m, CodifyVariables(o))
        @test association(est, x, y, z) isa Real
    end
end;