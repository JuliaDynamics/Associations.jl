using Random
rng = MersenneTwister(1234)

k = 5
ests = [
    GaussianMI(),
    KSG1(; k),
    KSG2(; k),
    GaoKannanOhViswanath(; k),
    GaoOhViswanath(; k),
]


x = StateSpaceSet(rand(rng, 1000, 1))
y = StateSpaceSet(rand(rng, 1000, 1))
@testset "MIShannon" begin
    @test MIShannon(base = 2) isa MIShannon

    # ----------------------------------------------------------------
    # Dedicated estimators.
    # ----------------------------------------------------------------
    # Just test that each estimator is reasonably close to zero for data from a uniform
    # distribution. This number varies wildly between estimators, so we're satisfied
    # to test just that they don't blow up.
    @testset "$(typeof(ests[i]).name.name)" for i in eachindex(ests)
        est = ests[i]
        measure = MIShannon(base = 2)
        mi = mutualinfo(measure, ests[i], x, y)
        @test mi isa Real
        @test -0.5 < mi < 0.1



        N = 5000
        c = 0.9
        Σ = [1 c; c 1]
        N2 = MvNormal([0, 0], Σ)
        mitrue = -0.5*log(det(Σ)) # in nats
        D2 = StateSpaceSet([rand(rng, N2) for i = 1:N])
        X = D2[:, 1] |> StateSpaceSet
        Y = D2[:, 2] |> StateSpaceSet

        mitrue_nats = -0.5*log(det(Σ))
        mitrue_bits = mitrue_nats / log(ℯ, 2)
        estimated_nats = mutualinfo(MIShannon(; base = ℯ), est, X, Y)
        estimated_bits = mutualinfo(MIShannon(; base = 2), est, X, Y)
        estimated_bits_kr = mutualinfo(MIShannon(; base = 2), Kraskov(), X, Y)

        @show mitrue_nats, mitrue_bits
        @show estimated_nats, estimated_bits, estimated_bits_kr
    end

end
