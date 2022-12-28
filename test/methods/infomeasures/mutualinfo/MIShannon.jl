using CausalityTools
using StateSpaceSets: Dataset

probests = [
    ValueHistogram(RectangularBinning(3))
    ValueHistogram(FixedRectangularBinning(0, 1, 3))
    NaiveKernel(0.2) # probably shouldn't be used.
]
xw
k = 5
diff_entropy_estimators = [
    Kraskov(; k),
    KozachenkoLeonenko(),
    GaoNaive(; k),
    GaoNaiveCorrected(; k),
    ZhuSingh(; k),
    Zhu(; k),
    Goria(; k),
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

    @testset "Definition: ShannonH3" begin
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
        # ----------------------------------------------------------------
        # We can't guarantee that the result is any particular value, because these are just
        # plug-in estimators. Just check that pluggin in works.
        @testset "$(typeof(probests[i]).name.name)" for i in eachindex(probests)
            def = ShannonH3()
            measure = MIShannon(base = 2)
            @test mutualinfo(def, measure, probests[i], x, y) isa Real
            @test mutualinfo(measure, probests[i], x, y) isa Real # default
        end

        # ----------------------------------------------------------------
        # Entropy-based estimators.
        # ----------------------------------------------------------------
        @testset "$(typeof(diff_entropy_estimators[i]).name.name)" for i in eachindex(diff_entropy_estimators)
            def = ShannonH3()
            measure = MIShannon(base = 2)
            mi = mutualinfo(measure, diff_entropy_estimators[i], x, y)
            @test mi isa Real
            @test -0.5 < mi < 0.1
        end
    end
end
