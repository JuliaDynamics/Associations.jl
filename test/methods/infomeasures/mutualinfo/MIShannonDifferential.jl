using CausalityTools
using StateSpaceSets: Dataset

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

# Test each definition separately.
@testset "MIShannonDifferential" begin
    @test MIShannonDifferential() isa MIShannonDifferential
    @test MIShannonDifferential(base = 2) isa MIShannonDifferential

    x = Dataset(rand(10000, 1))
    y = Dataset(rand(10000, 1))

    # Just test that each estimator is reasonably close to zero for data from a uniform
    # distribution. This number varies wildly between estimators, so we're satisfied
    # to test just that they don't blow up.
    @testset "$(typeof(diff_mi_estimators[i]).name.name)" for i in eachindex(diff_mi_estimators)
        measure = MIShannonDifferential(base = 2)
        est = diff_mi_estimators[i]
        mi = mutualinfo(measure, est, x, y)
        @test mi isa Real
        @test -0.5 < mi < 0.1
    end

    @testset "Definition: ShannonH3Differential" begin
        @testset "$(typeof(diff_entropy_estimators[i]).name.name)" for i in eachindex(diff_entropy_estimators)
            measure = MIShannonDifferential(base = 2, definition = ShannonH3Differential())
            est = diff_entropy_estimators[i]
            mi = mutualinfo(measure, est, x, y)
            @test mi isa Real
            @test -0.5 < mi < 0.1
        end
    end
end
