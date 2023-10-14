using Test
using CausalityTools
using Random
rng = MersenneTwister(1234)

# Double-sum estimation.
x = randn(rng, 50)
y = randn(rng, 50)
z = randn(rng, 50)

# Here we test all the possible "generic" ways of estimating `TEShannon`.
# Remaining tests are in the dedicated estimator test files, e.g. `Zhu1.jl`.
@test information(Zhu1(TEShannon(), k = 3), x, z) isa Real
@test information(Zhu1(TEShannon(), k = 3), x, z, y) isa Real
