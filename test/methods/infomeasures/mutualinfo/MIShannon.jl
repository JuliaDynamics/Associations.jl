using CausalityTools
using StateSpaceSets: Dataset

probests = [
    ValueHistogram(RectangularBinning(3))
    ValueHistogram(FixedRectangularBinning(0, 1, 3))
    NaiveKernel(0.2) # probably shouldn't be used.
]

@testset "MIShannon" begin
    @test MIShannon(base = 2) isa MIShannon

    x = Dataset(rand(10000, 2))
    y = Dataset(rand(10000, 1))

    @testset "Definition: ShannonH3" begin
        # We can't guarantee that the result is any particular value, because these are just
        # plug-in estimators. Just check that pluggin in works.
        @testset "$(typeof(probests[i]).name.name)" for i in eachindex(probests)
            measure = MIShannon(base = 2, definition = ShannonH3())
            mi = mutualinfo(measure, probests[i], x, y)
            @test mi isa Real
        end
    end
end
