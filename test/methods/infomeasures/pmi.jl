using StableRNGs
estd = Contingency(Dispersion(m = 2, c = 3))
esth = Contingency(ValueHistogram(5))
rng = StableRNG(123)

sys = system(Logistic4Chain(xi = rand(rng, 4); rng))
x, y, z, w = columns(trajectory(sys, 1000, Ttr = 10000))
@test estimate(PMI(), estd, x, w, z) >= 0
# Test that multivariate marginals work too.
@test estimate(PMI(), esth, x, w, Dataset(z, y)) >= 0

@test pmi(PMI(base = 3), esth, x, w, z) >= 0
@test pmi(SymbolicPermutation(m = 3), x, w, z) >= 0

x = rand(rng, 1:3, 20000)
y = rand(rng, 1:3, 20000)
z = rand(rng, 1:3, 20000)
@test round(pmi(CountOccurrences(), x, y, z), digits = 3) == 0
