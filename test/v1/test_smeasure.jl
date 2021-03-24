@testset "SMeasure" begin 
    x, y = rand(100), rand(100)
    X, Y = Dataset(rand(100, 3)), Dataset(rand(100, 2))
    dx, τx = 2, 1
    dy, τy = 2,1 

    @test s_measure(x, y) isa Float64
    @test s_measure(x, Y, dy = dy, τy = τy) isa Float64
    @test s_measure(X, y, dx = dx, τx = τx) isa Float64
    @test s_measure(X, Y) isa Float64
end