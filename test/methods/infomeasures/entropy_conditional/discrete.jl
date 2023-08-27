ests = [
    ValueHistogram(RectangularBinning(3)),
    ValueHistogram(FixedRectangularBinning(0:0.25:1.0)),
    OrdinalPatterns(m = 3),
    Dispersion(m = 2, c = 3)
]

# We don't have analytical tests for these estimators, so just test that
# they compute something.
x = rand(100)
y = rand(100)
@testset "$(typeof(ests[i]).name.name)" for i in eachindex(ests)
    est = ests[i]
    @test entropy_conditional(CEShannon(), est, x, y) >= 0
    @test entropy_conditional(CETsallisAbe(), est, x, y) isa Real
    @test_throws ArgumentError entropy_conditional(CETsallisFuruichi(), est, x, y)
end
