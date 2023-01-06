using CausalityTools
using StateSpaceSets: Dataset

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
    KSG1(; k),
    KSG2(; k),
    GaoKannanOhViswanath(; k),
    GaoOhViswanath(; k),
]


@testset "MIShannon" begin
    @test MIShannon(base = 2) isa MIShannon

    x = Dataset(rand(1000, 2))
    y = Dataset(rand(1000, 1))
    w1 = randn(10000)
    w2 = randn(10000)

    @test MIShannon(Shannon(; base = 2)) isa MIShannon
    @test MIShannon(base = 2) isa MIShannon

    # ----------------------------------------------------------------
    # Dedicated estimators.
    # ----------------------------------------------------------------
    # Just test that each estimator is reasonably close to zero for data from a uniform
    # distribution. This number varies wildly between estimators, so we're satisfied
    # to test just that they don't blow up.
    @testset "$(typeof(diff_mi_estimators[i]).name.name)" for i in eachindex(diff_mi_estimators)
        measure = MIShannon(base = 2)
        mi = mutualinfo(measure, diff_mi_estimators[i], x, y)
        @test mi isa Real
        @test -0.5 < mi < 0.1
    end

    # ----------------------------------------------------------------
    # Probability-based estimators.
    #
    # We can't guarantee that the result is any particular value, because these are just
    # plug-in estimators. Just check that pluggin in works.
    # ----------------------------------------------------------------

    # Estimators that accept dataset inputs
    @testset "$(typeof(probests[i]).name.name)" for i in eachindex(probests)
        m = MIShannon(base = 2)
        @test mutualinfo(m, probests[i], x, y) isa Real # default
    end

    @testset "$(typeof(probests_for_timeseries[i]).name)" for i in eachindex(probests_for_timeseries)
        est = probests_for_timeseries[i]
        measure = MIShannon(base = 2)
        @test mutualinfo(measure, est, w1, w2) isa Real # default
        @test mutualinfo(measure, est, w1, w2) >= 0 # default
        @show est
        # Doesn't work for datasets
        @test_throws MethodError mutualinfo(measure, est, x, y)
    end

    # ----------------------------------------------------------------
    # Entropy-based estimators.
    # ----------------------------------------------------------------
    @testset "$(typeof(diff_entropy_estimators[i]).name.name)" for i in eachindex(diff_entropy_estimators)
        measure = MIShannon(base = 2)
        mi = mutualinfo(measure, diff_entropy_estimators[i], x, y)
        @test mi isa Real
        @test -0.5 < mi < 0.15
    end
end
