using Random
rng = MersenneTwister(123456)

x, y = rand(rng, 100), rand(rng, 100)
@test crossmap(x, y, 3, 1) isa Float64
@test crossmap(x, y, 3, 1, :random) isa Vector{Float64}
@test crossmap(x, y, 3, 1, :segment) isa Vector{Float64}
