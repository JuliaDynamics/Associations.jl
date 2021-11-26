@testset "SMeasure" begin 
    x, y = rand(100), rand(100)
    X, Y = Dataset(rand(100, 3)), Dataset(rand(100, 2))
    Z, W = Dataset(rand(110, 2)), Dataset(rand(90, 4))
    dx, τx = 2, 1
    dy, τy = 2, 1 

    @test s_measure(x, y) isa Float64
    @test s_measure(x, Y, dx = dx, τx = τx) isa Float64
    @test s_measure(X, y, dy = dy, τy = τy) isa Float64
    @test s_measure(X, Y) isa Float64
    # test that multivariate datasets are being length-matched
    @test s_measure(X, Z) isa Float64
    @test s_measure(W, X) isa Float64
end