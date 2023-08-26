# Pre-discretized
x = rand(["a", "b", "c"], 200)
y = rand(["hello", "yoyo", "heyhey"], 200)
z = rand(["a", "b"], 200)
@test condmutualinfo(CMIShannon(), Contingency(), x, y, z) >= 0.0
@test condmutualinfo(CMIShannon(), Contingency(), x, y, z) isa Real
@test condmutualinfo(CMIRenyiJizba(), Contingency(), x, y, z) isa Real
@test condmutualinfo(CMIRenyiSarbu(), Contingency(), x, y, z) isa Real

# With discretization using a probabilities estimator
a, b, c = rand(100), rand(100), rand(100)
est = OrdinalPatterns(m = 2)
s1 = condmutualinfo(CMIShannon(), Contingency(est), a, b, c)
s2 = condmutualinfo(CMIShannon(), est, a, b, c)
@test s1 >= 0.0
@test s1 ≈ s2
@test condmutualinfo(CMIShannon(), Contingency(est), a, b, c) isa Real
@test condmutualinfo(CMIRenyiJizba(), Contingency(est), a, b, c) isa Real
@test condmutualinfo(CMIRenyiSarbu(), Contingency(est), a, b, c) isa Real
