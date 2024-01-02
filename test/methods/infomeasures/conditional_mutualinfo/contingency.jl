using Test
using Random
rng = Xoshiro(1234)

# Pre-discretized
x = rand(rng, ["a", "b", "c"], 200)
y = rand(rng, ["hello", "yoyo", "heyhey"], 200)
z = rand(rng, ["a", "b"], 200)
@test condmutualinfo(CMIShannon(), Contingency(), x, y, z) >= 0.0
@test condmutualinfo(CMIShannon(), Contingency(), x, y, z) isa Real
@test condmutualinfo(CMIRenyiJizba(), Contingency(), x, y, z) isa Real
@test condmutualinfo(CMIRenyiSarbu(), Contingency(), x, y, z) isa Real

# With discretization using a probabilities estimator
a, b, c = rand(rng, 100), rand(rng, 100), rand(rng, 100)
est = OrdinalPatterns{2}()
s1 = condmutualinfo(CMIShannon(), Contingency(est), a, b, c)
s2 = condmutualinfo(CMIShannon(), est, a, b, c)
@test s1 >= 0.0
@test s1 ≈ s2
@test condmutualinfo(CMIShannon(), Contingency(est), a, b, c) isa Real
@test condmutualinfo(CMIRenyiJizba(), Contingency(est), a, b, c) isa Real
@test condmutualinfo(CMIRenyiSarbu(), Contingency(est), a, b, c) isa Real
