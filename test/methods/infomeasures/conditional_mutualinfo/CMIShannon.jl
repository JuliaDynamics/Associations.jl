using CausalityTools
using StateSpaceSets: StateSpaceSet

probests = [
    ValueHistogram(RectangularBinning(3))
    #ValueHistogram(FixedRectangularBinning(0, 1, 3))
    NaiveKernel(0.2) # probably shouldn't be used.
]

probests_for_timeseries = [
    SymbolicPermutation(m = 3),
    Dispersion(c = 3, m = 2)
]

k = 5
diff_entropy_estimators = [
    Kraskov(; k),
    KozachenkoLeonenko(),
    ZhuSingh(; k),
    Zhu(; k),
    LeonenkoProzantoSavani(; k),
    Lord(; k = k*5),
]

diff_mi_estimators = [
    GaussianMI(),
    KSG1(; k),
    KSG2(; k),
    GaoKannanOhViswanath(; k),
    GaoOhViswanath(; k),
]


x = StateSpaceSet(rand(1000, 2))
y = StateSpaceSet(rand(1000, 1))
z = StateSpaceSet(rand(1000, 1))

@test FPVP() isa FPVP
@test MesnerShalisi() isa MesnerShalizi
@test MesnerShalizi() isa MesnerShalizi
@test PoczosSchneiderCMI() isa PoczosSchneiderCMI
@test Rahimzamani() isa Rahimzamani
@test GaussianCMI() isa GaussianCMI
@test condmutualinfo(FPVP(), x, y, z) isa Real
@test condmutualinfo(MesnerShalizi(), x, y, z) isa Real
@test condmutualinfo(PoczosSchneiderCMI(), x, y, z) isa Real
@test condmutualinfo(Rahimzamani(), x, y, z) isa Real
@test condmutualinfo(GaussianCMI(), x, y, z) isa Real

@test_throws ArgumentError condmutualinfo(CMIShannon(), FPVP(), x, y)
@test_throws ArgumentError condmutualinfo(CMIShannon(), FPVP(), x)

@testset "CMIShannon" begin
    @test m = CMIShannon(base = 2) isa CMIShannon

    x = StateSpaceSet(rand(10000, 2))
    y = StateSpaceSet(rand(10000, 1))
    z = StateSpaceSet(rand(10000, 1))
    w = StateSpaceSet(rand(10000, 1))

    @testset "Defaults" begin
        s, t, c = rand(100), rand(100), rand(100)
        est_bin = ValueHistogram(RectangularBinning(3))
        est_ksg = KSG1()

        # binning estimator yields non-negative values
        @test condmutualinfo(CMIShannon(base = 2), est_bin, s, t, c) >= 0.0
        @test condmutualinfo(CMIShannon(base = 2), est_ksg, s, t, c) isa Real # not guaranteed to be >= 0
    end

    @testset "Definition: CMIDefinitionShannonH4" begin
        @test CMIShannon() isa CMIShannon
        # ----------------------------------------------------------------
        # Dedicated estimators.
        # ----------------------------------------------------------------
        # Just test that each estimator is reasonably close to zero for data from a uniform
        # distribution. This number varies wildly between estimators, so we're satisfied
        # to test just that they don't blow up.
        @testset "$(typeof(diff_mi_estimators[i]).name.name)" for i in eachindex(diff_mi_estimators)
            est = diff_mi_estimators[i]
            mi = condmutualinfo(CMIShannon(base = 2), est, x, y, z)
            @test mi isa Real
            @test -0.5 < mi < 0.1
        end

        # ----------------------------------------------------------------
        # Probability-based estimators.
        #
        # We can't guarantee that the result is any particular value, because these are just
        # plug-in estimators. Just check that pluggin in works.
        # ----------------------------------------------------------------

        # Estimators that accept StateSpaceSet inputs
        @testset "$(typeof(probests[i]).name.name)" for i in eachindex(probests)
            est = probests[i]
            @test condmutualinfo(CMIShannon(base = 2), est, x, y, z) isa Real # default
        end

        # Estimators that only accept timeseries input
        a, b, c = rand(10000), rand(10000), rand(10000)

        @testset "$(typeof(probests_for_timeseries[i]).name)" for i in eachindex(probests_for_timeseries)
            est = probests_for_timeseries[i]
            cmi = CMIShannon(base = 2)
            @test condmutualinfo(cmi, est, a, b, c) >= 0
            @test condmutualinfo(cmi, est, x, y, z) >= 0
        end

        # ----------------------------------------------------------------
        # Entropy-based estimators.
        # ----------------------------------------------------------------
        @testset "$(typeof(diff_entropy_estimators[i]).name.name)" for i in eachindex(diff_entropy_estimators)
            est = diff_entropy_estimators[i]
            mi = condmutualinfo(CMIShannon(base = 2), est, x, y, z)
            @test mi isa Real
            @test -0.5 < mi < 0.1
        end
    end
end
