using CausalityTools
using StateSpaceSets: StateSpaceSet
using Test
using Random; rng = Xoshiro(1235)


x = StateSpaceSet(rand(rng, 1000, 2))
y = StateSpaceSet(rand(rng, 1000, 1))
z = StateSpaceSet(rand(rng, 1000, 1))

@testset "Dedicated estimators" begin
    @test FPVP() isa FPVP
    @test MesnerShalisi() isa MesnerShalizi
    @test MesnerShalizi() isa MesnerShalizi
    @test PoczosSchneiderCMI() isa PoczosSchneiderCMI
    @test Rahimzamani() isa Rahimzamani
    @test GaussianCMI() isa GaussianCMI

    @test association(FPVP(), x, y, z) isa Real
    @test association(MesnerShalizi(), x, y, z) isa Real
    @test association(PoczosSchneiderCMI(), x, y, z) isa Real
    @test association(Rahimzamani(), x, y, z) isa Real
    @test association(GaussianCMI(), x, y, z) isa Real
end

@testset "Input checks" begin
    @test_throws ArgumentError association(CMIShannon(), FPVP(), x, y)
    @test_throws ArgumentError association(CMIShannon(), FPVP(), x)
end

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

@testset "EntropyDecomposition" begin

    @testset "Discrete" begin
        outcome_spaces = [
            ValueBinning(RectangularBinning(3)),
            OrdinalPatterns(m = 3),
            Dispersion(c = 3, m = 2)
        ]

        for (i, o) in enumerate(outcome_spaces)
            est = EntropyDecomposition(est_diff, o)
            @test association(est, x, y, z) ≥ 0
        end
    end

    @testset "Differential" begin
        
        diff_entropy_estimators = [
            Kraskov(; k),
            KozachenkoLeonenko(),
            ZhuSingh(; k),
            Zhu(; k),
            LeonenkoProzantoSavani(; k),
            Lord(; k = k*5),
        ]
        for (i, est_diff) in enumerate(diff_entropy_estimators)
            @test association(est_diff, x, y, z) isa Real # not guaranteed to be >= 0
        end
    end
end

@testset "MIDecomposition" begin
    k = 2
    diff_mi_estimators = [
        GaussianMI(),
        KSG1(; k),
        KSG2(; k),
        GaoKannanOhViswanath(; k),
        GaoOhViswanath(; k),
    ]

    for (i, est_mi) in enumerate(diff_mi_estimators)
        est = MIDecomposition(CMIShannon(), est_mi)
        @test association(est, x, y, z) isa Real # not guaranteed to be >= 0
    end
end