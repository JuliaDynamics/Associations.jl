
@testset "Leanings" begin
    x = [0, 0, 1, 0, 0, 1, 0, 0, 1, 0]
    y = [0, 0, 0, 1, 0, 0, 1, 0, 0, 1]
    @test penchant(x, y, 1, weighted = false) ≈ 1.0
    @test penchant(y, x, 1, weighted = false) ≈ 1/7

    @test penchant(x, y, 1, weighted = true) ≈ 1.0
    @test penchant(y, x, 1, weighted = true) ≈ 3/63

    @test lean(x, y, 1, weighted = true) ≈ 60/63
end