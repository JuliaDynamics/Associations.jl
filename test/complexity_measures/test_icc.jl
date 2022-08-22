@testset "ICC" begin
    x = rand(0:1, 500)
    y = rand(0:1, 500)
    est = CompressionComplexityCausalityEstimator(
        algorithm = EffortToCompress(normalize = true),
        w = 15,
        L = 30,
        step = 10)
    @test icc(x, y, est) isa typeof(1.0)
end
