using StableRNGs
estd = Contingency(Dispersion(m = 2, c = 3))
esth = Contingency(ValueHistogram(5))
rng = StableRNG(123)

sys = system(Logistic4Chain(xi = rand(rng, 4); rng))
x, y, z, w = columns(trajectory(sys, 1000, Ttr = 10000))
@test estimate(PMI(), estd, x, w, z) >= 0
# Test that multivariate marginals work too.
@test estimate(PMI(), esth, x, w, Dataset(z, y)) >= 0
