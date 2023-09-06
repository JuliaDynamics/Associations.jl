ests = [
    ValueBinning(RectangularBinning(3)),
    ValueBinning(FixedRectangularBinning(0:0.25:1.0)),
    OrdinalPatterns(m = 3),
    Dispersion(m = 2, c = 3)
]

# We don't have analytical tests for these estimators, so just test that
# they compute something.
x = rand(100)
y = rand(100)
@testset "$(typeof(ests[i]).name)" for i in eachindex(ests)
    est = ests[i]
    @test mutualinfo(MIShannon(), est, x, y) >= 0
    @test mutualinfo(MITsallisFuruichi(), est, x, y) isa Real
    @test mutualinfo(MITsallisMartin(), est, x, y) isa Real
    @test mutualinfo(MIRenyiJizba(), est, x, y) isa Real
    @test_throws ArgumentError mutualinfo(MIRenyiSarbu(), est, x, y) isa Real
end
