using Random
rng = MersenneTwister(123456)

x, y = rand(rng, 1000), rand(rng, 1000)
@test pai(x, y, 3, 1) isa Float64
@test pai(x, y, 3, 1, :random) isa Vector{Float64}
@test pai(x, y, 3, 1, :segment) isa Vector{Float64}
