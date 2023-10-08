using Test
using CausalityTools
using Random
rng = MersenneTwister(1234)

# Double-sum estimation.
x = rand(rng, ["a", "b", "c"], 200)
y = rand(rng, ["hello", "yoyo", "heyhey"], 200)

# The estimation of probabilities is decoupled from the estimation of the mutual info.
# We could in principle use any probabilities estimator here, but we default to `RelativeAmount`.
p = probabilities(x, y) 
@test information(MIShannon(), p) >= 0.0

# TODO: three-entropies sum