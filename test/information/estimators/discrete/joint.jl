using CausalityTools
using Test
using Random
rng = Xoshiro(1234)

# ------------------------------------------------------------
# API test 
# ------------------------------------------------------------
x = randn(rng, 1000)
y = x .+ randn(rng, 1000)
z = y .+ randn(rng, 1000)

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

@testset "JointProbabilities estimator with $m" for m in bivariate_symmetric_measures
    @testset "CodifyVariables with $o" for o in os
        est = JointProbabilities(m, CodifyVariables(o))
        est_xy = information(est, x, y)
        est_yx = information(est, y, x)

        @test est_xy isa Real
        @test est_yx isa Real
        @test est_xy == est_yx # symmetry
    end
end

bivariate_asymmetric_measures = [
    CEShannon(),
    CETsallisAbe(),
    CETsallisFuruichi(),
    HellingerDistance(),
    KLDivergence(),
    RenyiDivergence(),
    VariationDistance(),
]

@testset "JointProbabilities estimator with $m" for m in bivariate_asymmetric_measures
    @testset "CodifyVariables with $o" for o in os
        est = JointProbabilities(m, CodifyVariables(o))
        @test information(est, x, y) isa Real
    end
end

trivariate_asymmetric_measures = [
    CMIShannon(),
    CMITsallis(),
    CMIRenyiSarbu(),
]

@testset "JointProbabilities estimator with $m" for m in trivariate_asymmetric_measures
    @testset "CodifyVariables with $o" for o in os
        est = JointProbabilities(m, CodifyVariables(o))
        @test information(est, x, y, z) isa Real
    end
end