using Entropies: Kraskov, KozachenkoLeonenko, VisitationFrequency, RectangularBinning
@testset "Mutual information" begin
    x = rand(100)
    y = rand(100)
    z = Dataset(rand(100, 2))
    w = Dataset(rand(100, 3))

    # Estimators for which Renyi entropies can be used
    est = VisitationFrequency(RectangularBinning(0.2))
    @test mutualinfo(est, x, y) isa Real
    @test mutualinfo(est, x, y) isa Real
    @test mutualinfo(est, x, z) isa Real
    @test mutualinfo(est, z, x) isa Real
    @test mutualinfo(est, z, w) isa Real

    # Estimators for which Renyi entropies cannot be used
    est_kl = KozachenkoLeonenko()
    est_k = Kraskov(k = 2)
    est_k1 = KSG1(k = 2)
    est_k2 = KSG2(k = 2)

    @test mutualinfo(est_kl, x, y) isa Real
    @test mutualinfo(est_kl, x, z) isa Real
    @test mutualinfo(est_kl, z, x) isa Real
    @test mutualinfo(est_kl, z, w) isa Real

    @test mutualinfo(est_k, x, y) isa Real
    @test mutualinfo(est_k, x, z) isa Real
    @test mutualinfo(est_k, z, x) isa Real
    @test mutualinfo(est_k, z, w) isa Real

    @test mutualinfo(est_k1, x, y) isa Real
    @test mutualinfo(est_k1, x, z) isa Real
    @test mutualinfo(est_k1, z, x) isa Real
    @test mutualinfo(est_k1, z, w) isa Real

    @test mutualinfo(est_k2, x, y) isa Real
    @test mutualinfo(est_k2, x, z) isa Real
    @test mutualinfo(est_k2, z, x) isa Real
    @test mutualinfo(est_k2, z, w) isa Real
end
