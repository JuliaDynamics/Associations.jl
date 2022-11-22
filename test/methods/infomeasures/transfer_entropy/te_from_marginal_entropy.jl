using Entropies: Zhu, ZhuSingh, Kraskov, KozachenkoLeonenko
using Entropies: ValueHistogram, TransferOperator, RectangularBinning
using Random
rng = MersenneTwister(123456)

n = 50000
# TE(x -> y) and TE(x -> y | z) should be zero for independent data.
x, y, z = randn(rng, n), randn(rng, n), randn(rng, n)
base = 2
@testset "Marginal estimation (indirect entropies)" begin
    indirect_estimators = [
        Zhu(; k = 5, base),
        ZhuSingh(; k = 5, base),
        Kraskov(; k = 5, base),
        KozachenkoLeonenko(; base)
    ]

    @testset "$(indirect_estimators[i])" for i in eachindex(indirect_estimators)
        e = indirect_estimators[i]
        te = transferentropy(e, x, y)
        tec = transferentropy(e, x, y, z)
        # Check that we're "sufficiently close" to zero. This is completely heuristic.
        @test round(abs(te), digits = 2) <= 0.03
        @test round(abs(tec), digits = 2) <= 0.03
    end
end

@testset "Marginal estimation (direct entropies)" begin
    direct_estimators = [
        ValueHistogram(RectangularBinning(3)),
        TransferOperator(RectangularBinning(3)),
        NaiveKernel(0.5),
    ]
    e = Shannon(; base)
    @testset "$(direct_estimators[i])" for i in eachindex(direct_estimators)
        est = direct_estimators[i]
        te = transferentropy(e, est, x, y)
        tec = transferentropy(e, est, x, y, z)
        # Check that we're "sufficiently close" to zero. This is completely heuristic.
        # Different estimators have different bias, so we need to individually adjust
        # bounds here.
        @test round(abs(te), digits = 2) <= 0.03
        @test round(abs(tec), digits = 2) <= 0.03
    end
end
