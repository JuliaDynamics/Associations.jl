using Test
using CausalityTools

# There should be zero information gain from `x` over `y` for independent random variables.
using Random
rng = Xoshiro(1234)
n = 100000
x, y = rand(rng, n), rand(rng, n)
d = CodifyVariables(OrdinalPatterns(m = 3))
div_kl = association(JointProbabilities(KLDivergence(), d), x, y)
@test abs(div_kl) ≤ 0.001
