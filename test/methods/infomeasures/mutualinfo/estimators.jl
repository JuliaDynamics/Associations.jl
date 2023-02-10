ests = [
    GaussianMI(),
    KSG1(; k),
    KSG2(; k),
    GaoKannanOhViswanath(; k),
    GaoOhViswanath(; k),
]


x = Dataset(rand(300, 2))
y = Dataset(rand(300, 1))
@testset "MIShannon" begin
    @test MIShannon(base = 2) isa MIShannon

    # ----------------------------------------------------------------
    # Dedicated estimators.
    # ----------------------------------------------------------------
    # Just test that each estimator is reasonably close to zero for data from a uniform
    # distribution. This number varies wildly between estimators, so we're satisfied
    # to test just that they don't blow up.
    @testset "$(typeof(ests[i]).name.name)" for i in eachindex(ests)
        measure = MIShannon(base = 2)
        mi = mutualinfo(measure, ests[i], x, y)
        @test mi isa Real
        @test -0.5 < mi < 0.1
    end

end
