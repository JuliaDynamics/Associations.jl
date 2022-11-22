using Entropies: VisitationFrequency, RectangularBinning
using DelayEmbeddings: Dataset

@testset "Conditional mutual information" begin
    s, t, c = rand(100), rand(100), rand(100)
    est_bin = ValueHistogram(RectangularBinning(3))
    # binning estimator yields non-negative values
    @test conditional_mutualinfo(est_bin, s, t, c) isa Real
    @test conditional_mutualinfo(est_bin, s, t, c) >= 0.0
    # verify formula I(X, Y | Z) = I(X; Y, Z) - I(X, Z)
    @test conditional_mutualinfo(est_bin, s, t, c) ≈
        mutualinfo(est_bin, s, Dataset(t, c)) - mutualinfo(est_bin, s, c)

    @test conditional_mutualinfo(Kraskov1(k = 2), s, t, c) isa Real
    @test conditional_mutualinfo(Kraskov1(k = 2), s, t, c) ≈
        mutualinfo(Kraskov1(k = 2), s, Dataset(t, c)) -
            mutualinfo(Kraskov1(k = 2), s, c)

    # Different types of input
    @test conditional_mutualinfo(est_bin, s, Dataset(t, c), c) isa Real
    @test conditional_mutualinfo(est_bin, Dataset(s, t), Dataset(t, c), c) isa Real
    @test conditional_mutualinfo(est_bin, Dataset(s, t), Dataset(t, c), Dataset(c, s)) isa Real
    @test conditional_mutualinfo(est_bin, s, Dataset(t, c), Dataset(c, s)) isa Real
    @test conditional_mutualinfo(est_bin, s, t, Dataset(c, s)) isa Real
    @test conditional_mutualinfo(est_bin, Dataset(s, t), t, c) isa Real
end
