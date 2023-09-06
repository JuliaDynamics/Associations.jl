using Test
using CausalityTools

# There should be zero information gain from `x` over `y` for independent random variables.
using Random
rng = Xoshiro(1234)
n = 100000
x, y = rand(rng, n), rand(rng, n)
div_kl = information(KLDivergence(), OrdinalPatterns(m=3), x, y)
@test abs(div_kl) ≤ 0.001
