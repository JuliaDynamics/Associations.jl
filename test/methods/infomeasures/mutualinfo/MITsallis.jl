using CausalityTools
using StateSpaceSets: Dataset

probests = [
    ValueHistogram(RectangularBinning(3))
    ValueHistogram(FixedRectangularBinning(0, 1, 3))
    NaiveKernel(0.2) # probably shouldn't be used.
]

probests_for_timeseries = [
    SymbolicPermutation(m = 3),
    Dispersion(c = 3, m = 2)
]

@testset "MITsallis" begin
    @test MITsallis(base = 2) isa MITsallis

    x = Dataset(rand(10000, 2))
    y = Dataset(rand(10000, 1))

    w1 = rand(10000)
    w2 = rand(10000)

    @testset "MITsallisFuruichi" begin
        # We can't guarantee that the result is any particular value, because these are just
        # plug-in estimators. Just check that pluggin in works.
        @testset "$(typeof(probests[i]).name.name)" for i in eachindex(probests)
            est = probests[i]
            mi_furu = MITsallisFuruichi(q = 1.5, base = 2)
            @test mutualinfo(mi_furu, est, x, y) isa Real
        end

        # Estimators that only accept timeseries input
        @testset "$(typeof(probests_for_timeseries[i]).name)" for i in eachindex(probests_for_timeseries)
            est = probests_for_timeseries[i]
            mi_furu = MITsallisFuruichi(q = 1.5, base = 2)
            @test mutualinfo(mi_furu, est, w1, w2) isa Real
            # Doesn't work for datasets.
            #@test_throws MethodError mutualinfo(measure, est, x, y)
        end
    end

    @testset "MITsallisMartin" begin
        # We can't guarantee that the result is any particular value, because these are just
        # plug-in estimators. Just check that pluggin in works.
        @testset "$(typeof(probests[i]).name.name)" for i in eachindex(probests)
            est = probests[i]
            mi_mart = MITsallisMartin(q = 1.5, base = 2)
            @test mutualinfo(mi_mart, est, x, y) isa Real
        end

        # Estimators that only accept timeseries input
        @testset "$(typeof(probests_for_timeseries[i]).name)" for i in eachindex(probests_for_timeseries)
            est = probests_for_timeseries[i]
            mi_mart = MITsallisMartin(q = 1.5, base = 2)
            @test mutualinfo(mi_mart, est, w1, w2) isa Real
            # Doesn't work for datasets.
            #@test_throws MethodError mutualinfo(measure, est, x, y)
        end
    end
end
