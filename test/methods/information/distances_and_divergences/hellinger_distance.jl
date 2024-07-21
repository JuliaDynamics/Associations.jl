using Test
using CausalityTools

# There should be zero information gain from `x` over `y` for independent random variables.
using Random
rng = Xoshiro(1234)
n = 100000
x, y = rand(rng, n), rand(rng, n)
disc = CodifyVariables(OrdinalPatterns(m = 3))
hel = association(JointProbabilities(HellingerDistance(), disc), x, y)
@test abs(hel) ≤ 0.001
