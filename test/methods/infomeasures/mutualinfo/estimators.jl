using Random
import Random: seed!
rng = MersenneTwister(1234)
using Distributions: MvNormal
using LinearAlgebra: det

k = 5
ests_mi = [
    GaussianMI(normalize=true),
    GaussianMI(normalize=false),
    KSG1(; k),
    KSG2(; k),
    GaoKannanOhViswanath(; k),
    GaoOhViswanath(; k),
]

ests_diffent = [
    Kraskov(; k),
    KozachenkoLeonenko(;),
    Zhu(; k),
    ZhuSingh(; k),
    Gao(; k),
    Goria(; k),
    Lord(; k)
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
    @testset "$(typeof(ests_mi[i]).name.name)" for i in eachindex(ests_mi)
        est = ests_mi[i]
        measure = MIShannon(base = 2)
        mi = mutualinfo(measure, est, x, y)
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
        mitrue_bits = CausalityTools._convert_logunit(mitrue_nats, ℯ, 2)
        estimated_nats = mutualinfo(MIShannon(; base = ℯ), est, X, Y)
        estimated_bits = mutualinfo(MIShannon(; base = 2), est, X, Y)
        estimated_bits_kr = mutualinfo(MIShannon(; base = 2), Kraskov(), X, Y)
    end

    # ----------------------------------------------------------------
    # `DifferentialInformationEstimator`s`
    # ----------------------------------------------------------------
    @testset "$(typeof(ests_diffent[i]).name.name)" for i in eachindex(ests_diffent)
        est = ests_diffent[i]
        measure = MIShannon(base = 2)
        mi = mutualinfo(measure, est, x, y)
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
        mitrue_bits = CausalityTools._convert_logunit(mitrue_nats, ℯ, 2)
        estimated_nats = mutualinfo(MIShannon(; base = ℯ), est, X, Y)
        estimated_bits = mutualinfo(MIShannon(; base = 2), est, X, Y)
        estimated_bits_kr = mutualinfo(MIShannon(; base = 2), Kraskov(), X, Y)
    end

end

@testset "GaussianMI" begin
    # The other estimator tests only compute whether the estimators run "at all".
    # For some special cases of the Gaussian we can also compare with a closed form solution.

    @testset "Normalized equals unnormalized" begin
        x′ = StateSpaceSet(2. .* x.data .+ [SVector(1.)])
        y′ = StateSpaceSet(3. .* y.data .- [SVector(1.)])
        @test (  mutualinfo(GaussianMI(normalize=false), x , y)
                 ≈ mutualinfo(GaussianMI(normalize=true) , x′, y′))
    end

    @testset "Compare with analytic eq" begin
        # Test based on https://en.wikipedia.org/wiki/Mutual_information#Linear_correlation.
        # We choose parameters arbitrarily:
        @testset "normalize=false" begin
            σ_1 = 1.0
            σ_2 = 1.0
            ρ = 0.5

            μ = [0.; 0.]
            Σ = [σ_1^2 ρ*σ_1*σ_2;
                ρ*σ_1*σ_2 σ_2^2]

            seed!(rng, 1)
            xys = rand(rng, MvNormal(μ, Σ), 1_000_000)
            # Notice that MIShannon.base is `2` by default, but math expects `ℯ`.
            @test estimate(
                MIShannon(; base=ℯ),
                GaussianMI(normalize=false),
                StateSpaceSet(xys[1:1, :]'),
                StateSpaceSet(xys[2:2, :]')
            ) ≈ -1/2 * log(1 - ρ^2)  atol=1e-3
        end
        @testset "normalize=true" begin
            σ_1 = 0.5
            σ_2 = 1.5
            ρ = 0.5

            μ = [1.5; 2.5]
            Σ = [σ_1^2 ρ*σ_1*σ_2;
                ρ*σ_1*σ_2 σ_2^2]

            seed!(rng, 1)
            xys = rand(rng, MvNormal(μ, Σ), 1_000_000)
            # Notice that MIShannon.base is `2` by default, but math expects `ℯ`.
            @test estimate(
                MIShannon(; base=ℯ),
                GaussianMI(normalize=true),
                StateSpaceSet(xys[1:1, :]'),
                StateSpaceSet(xys[2:2, :]')
            ) ≈ -1/2 * log(1 - ρ^2)  atol=1e-3
        end
    end
end
