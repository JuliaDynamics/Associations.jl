using StableRNGs
rng = StableRNG(123)

ed = Dispersion(m = 2, c = 2)
vh = ValueBinning(2)
sp = OrdinalPatterns(m=2)
estd = Contingency(ed)
esth = Contingency(vh)
ests = Contingency(sp)

sys = system(Logistic4Chain(xi = rand(rng, 4); rng))
x, y, z, w = columns(first(trajectory(sys, 50, Ttr = 10000)))
ZW = StateSpaceSet(z, w)
@test pmi(estd, x, y, z) >= 0
@test pmi(esth, x, y, z) >= 0
@test pmi(ests, x, y, z) >= 0
@test pmi(estd, x, y, z) >= 0
@test pmi(ed, x, y, z) >= 0
@test pmi(vh, x, y, z) >= 0
@test pmi(sp, x, y, z) >= 0
@test pmi(estd, x, y, ZW) >= 0
@test pmi(esth, x, y, ZW) >= 0
@test pmi(ests, x, y, ZW) >= 0
@test pmi(estd, x, y, ZW) >= 0
@test pmi(ed, x, y, ZW) >= 0
@test pmi(vh, x, y, ZW) >= 0
@test pmi(sp, x, y, ZW) >= 0


sys = system(Logistic4Chain(xi = rand(rng, 4); rng))
x, y, z, w = columns(first(trajectory(sys, 1000, Ttr = 10000)))
@test estimate(PMI(), estd, x, w, z) >= 0
# Test that multivariate marginals work too.
@test estimate(PMI(), esth, x, w, Dataset(z, y)) >= 0

@test pmi(PMI(base = 3), esth, x, w, z) >= 0
@test pmi(OrdinalPatterns{3}(), x, w, z) >= 0

x = rand(rng, 1:3, 20000)
y = rand(rng, 1:3, 20000)
z = rand(rng, 1:3, 20000)
@test round(pmi(UniqueElements(), x, y, z), digits = 3) == 0

# Independence tests
x = rand(rng, 50)
y = rand(rng, 50)
z = rand(rng, 50)
X = StateSpaceSet(x)
Y = StateSpaceSet(y)
Z = StateSpaceSet(z)

nshuffles = 5
lptest = LocalPermutationTest(PMI(), OrdinalPatterns(); nshuffles, rng)
sutest = SurrogateTest(PMI(), OrdinalPatterns(); nshuffles, rng)
@test independence(lptest, x, y, z) isa LocalPermutationTestResult
@test independence(lptest, X, Y, Z) isa LocalPermutationTestResult
@test independence(sutest, x, y, z) isa SurrogateTestResult
@test independence(sutest, X, Y, Z) isa SurrogateTestResult
