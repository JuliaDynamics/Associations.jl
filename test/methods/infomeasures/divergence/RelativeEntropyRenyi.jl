using CausalityTools
using StateSpaceSets: Dataset

probests = [
    ValueHistogram(RectangularBinning(3))
    ValueHistogram(FixedRectangularBinning(0, 1, 3))
]

@testset "RelativeEntropyRenyi" begin
    @test RelativeEntropyRenyi(base = 2) isa RelativeEntropyRenyi

    x = Dataset(rand(10000, 2))
    y = Dataset(rand(10000, 2))

    @testset "Definition: RenyiDivergence" begin
        # We can't guarantee that the result is any particular value, because these are just
        # plug-in estimators. Just check that pluggin in works.
        @testset "$(typeof(probests[i]).name.name)" for i in eachindex(probests)
            measure = RelativeEntropyRenyi(base = 2, definition = RenyiDivergence())
            est = probests[i]
            dv = divergence(measure, est, x, y)
            @test dv isa Real
        end
    end
end
