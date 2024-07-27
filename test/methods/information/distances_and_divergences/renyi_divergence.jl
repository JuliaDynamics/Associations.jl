using Test
using Associations

# There should be zero information gain from `x` over `y` for independent random variables.
using Random
rng = Xoshiro(1234)
n = 100000
x, y = rand(rng, n), rand(rng, n)
d = CodifyVariables(OrdinalPatterns(m = 3))

@test div_r = association(JointProbabilities(RenyiDivergence(q = 0.5), d), x, y) ≤ 0.001
@test div_r = association(JointProbabilities(RenyiDivergence(q = 2), d), x, y) ≤ 0.001
@test div_r = association(JointProbabilities(RenyiDivergence(q = 4), d), x, y) ≤ 0.001
@test association(JointProbabilities(RenyiDivergence(q = Inf), d), x, y) ≥ 0.0
