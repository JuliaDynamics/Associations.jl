k = 5
ests = [
    Kraskov(; k),
    KozachenkoLeonenko(),
    ZhuSingh(; k),
    Zhu(; k),
    Lord(; k = k*5),
]

# We don't have analytical tests for these estimators, so just test that
# they compute something.
x = rand(100)
y = rand(100)
@testset "$(typeof(ests[i]).name)" for i in eachindex(ests)
    est = ests[i]
    @test association(MIShannon(), est, x, y) isa Real
    @test_throws ArgumentError association(MITsallisFuruichi(), est, x, y) isa Real
    @test_throws ArgumentError association(MITsallisMartin(), est, x, y) isa Real
    @test_throws ArgumentError association(MIRenyiJizba(), est, x, y) isa Real
    @test_throws ArgumentError association(MIRenyiSarbu(), est, x, y) isa Real
end
