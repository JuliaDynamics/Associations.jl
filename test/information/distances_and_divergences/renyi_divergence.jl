using Test
using CausalityTools

# There should be zero information gain from `x` over `y` for independent random variables.
using Random
rng = Xoshiro(1234)
n = 100000
x, y = rand(rng, n), rand(rng, n)
o = OrdinalPatterns(m=3)
@test div_r = information(RenyiDivergence(q = 0.5), o, x, y) ≤ 0.001
@test div_r = information(RenyiDivergence(q = 2), o, x, y) ≤ 0.001
@test div_r = information(RenyiDivergence(q = 4), o, x, y) ≤ 0.001
@test information(RenyiDivergence(q = Inf), o, x, y) ≥ 0.0
