using CausalityTools
using StateSpaceSets: Dataset

probests = [
    ValueHistogram(RectangularBinning(3))
    ValueHistogram(FixedRectangularBinning(0, 1, 3))
    NaiveKernel(0.2) # probably shouldn't be used.
]

@testset "MITsallis" begin
    @test MITsallis(base = 2) isa MITsallis

    x = Dataset(rand(10000, 2))
    y = Dataset(rand(10000, 1))

    @testset "Definition: TsallisH3" begin
        # We can't guarantee that the result is any particular value, because these are just
        # plug-in estimators. Just check that pluggin in works.
        @testset "$(typeof(probests[i]).name.name)" for i in eachindex(probests)
            measure = MITsallis(q = 1.5, base = 2, definition = TsallisH3())
            mi = mutualinfo(measure, probests[i], x, y)
            @test mi isa Real
        end
    end
end
