using CausalityTools
using StateSpaceSets: Dataset
using Distributions: MvNormal
using Random; rng = MersenneTwister(1234567)

probests = [
    ValueHistogram(RectangularBinning(3))
    ValueHistogram(FixedRectangularBinning(0, 1, 3))
]

@testset "RelativeEntropyRenyi" begin
    @test RelativeEntropyRenyi(base = 2) isa RelativeEntropyRenyi

    σx = 0.5; σy = 0.5; vx = 1.0; vy = 1.0
    Nx = MvNormal([0, 0], [vx σx; σx vx])
    Ny = MvNormal([0, 0], [vy σy; σy vy])
    x = Dataset([rand(rng, Nx) for i = 1:100000])
    y = Dataset([rand(rng, Ny) for i = 1:100000])

    @testset "Definition: RenyiDivergence" begin
        # We can't guarantee that the result is any particular value, because these are just
        # plug-in estimators. Just check that pluggin in works.
        @testset "$(typeof(probests[i]).name.name)" for i in eachindex(probests)
            def = RenyiDivergence()
            measure = RelativeEntropyRenyi(base = 2)
            est = probests[i]
            @test divergence(measure, est, x, y) isa Real
        end
    end
end
