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
    @test entropy_conditional(ConditionalEntropyShannon(), est, x, y) isa Real
    @test_throws ArgumentError  entropy_conditional(ConditionalEntropyTsallisAbe(), est, x, y)
    @test_throws ArgumentError entropy_conditional(ConditionalEntropyTsallisFuruichi(), est, x, y)
end
