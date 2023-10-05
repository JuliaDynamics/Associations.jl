using Statistics
@testset "Dynamical complexity" begin
    # Some random binary time series segments
    x, y = rand(0:1, 100), rand(0:1, 100)

    # Split segments into past and present
    X⁻, Xⁱ = x[1:49], x[50:end]
    alg = EffortToCompress(normalize = true)
    @test dynamical_complexity(Xⁱ, X⁻, alg) isa typeof(1.0)

    # Split segments into past and present
    alg = EffortToCompress(normalize = true)
    w, L, step = 10, 15, 5
    @test dynamical_complexity(Statistics.mean, x, y, alg, w, L, step) isa typeof(1.0)
end
