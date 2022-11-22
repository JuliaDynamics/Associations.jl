using Random
rng = Random.default_rng()
x = sin.(1:100) .+ rand(rng, 100)
y = sin.(1:100) .+ rand(rng, 100)
z = sin.(1:100) .+ rand(rng, 100)

@testset "Transfer entropy (Hilbert)" begin
    est_vf = ValueHistogram(RectangularBinning(3))
    special_estimators = [
        Hilbert(source = Phase(), target = Amplitude(), est_vf),
        Hilbert(source = Amplitude(), target = Amplitude(), est_vf),
        Hilbert(source = Phase(), target = Phase(), est_vf),
        Hilbert(source = Phase(), target = Phase(), cond = Amplitude(), est_vf),
    ]
    e = Shannon(; base)
    @testset "$(special_estimators[i])" for i in eachindex(special_estimators)
        est = special_estimators[i]
        te = transferentropy(e, est, x, y)
        tec = transferentropy(e, est, x, y, z)
        # Check that we're "sufficiently close" to zero. This is completely heuristic.
        # Different estimators have different bias, so we need to individually adjust
        # bounds here.
        @test round(abs(te), digits = 2) <= 0.03
        @test round(abs(tec), digits = 2) <= 0.03
    end
end
