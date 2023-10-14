using Test
using CausalityTools
using Random
rng = MersenneTwister(1234)

# Double-sum estimation.
x = randn(rng, 50)
y = randn(rng, 50)
z = randn(rng, 50)

# Here we test all the possible "generic" ways of estimating `TEShannon`.
# Remaining tests are in the dedicated estimator test files, e.g. `Zhu1.jl`.
est_diff = EntropyDecomposition(TEShannon(), Kraskov(k=3))
@test information(est_diff, x, z) isa Real
@test information(est_diff, x, z, y) isa Real

est_disc = EntropyDecomposition(TEShannon(), PlugIn(Shannon()), ValueBinning(2));
@test information(est_disc, x, z) isa Real
@test information(est_disc, x, z, y) isa Real

est_mi = MIDecomposition(TEShannon(), KSG1())
@test information(est_mi, x, z) isa Real
@test information(est_mi, x, z, y) isa Real

est_cmi = CMIDecomposition(TEShannon(), FPVP())
@test information(est_cmi, x, z) isa Real
@test information(est_cmi, x, z, y) isa Real
