# Pre-discretized
x = rand(["a", "b", "c"], 200)
y = rand(["hello", "yoyo", "heyhey"], 200)
@test mutualinfo(MIShannon(), Contingency(), x, y) >= 0.0

@test mutualinfo(MITsallisFuruichi(), Contingency(), x, y) isa Real
@test mutualinfo(MITsallisMartin(), Contingency(), x, y) isa Real
@test mutualinfo(MIRenyiJizba(), Contingency(), x, y) isa Real
@test mutualinfo(MIRenyiSarbu(), Contingency(), x, y) isa Real

# With discretization using a probabilities estimator
z, w = rand(100), rand(100)
est = OrdinalPatterns(m = 3)
@test mutualinfo(MIShannon(), Contingency(est), z, w) >= 0.0
@test mutualinfo(MITsallisFuruichi(), Contingency(est), z, w) isa Real
@test mutualinfo(MITsallisMartin(), Contingency(est), z, w) isa Real
@test mutualinfo(MIRenyiJizba(), Contingency(est), z, w) isa Real
