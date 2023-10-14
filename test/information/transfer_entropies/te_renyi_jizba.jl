using Test
using CausalityTools
using Random
rng = MersenneTwister(1234)

# Double-sum estimation.
x = randn(rng, 100)
y = randn(rng, 100)
z = randn(rng, 100)

def = TERenyiJizba(base = 3, q = 0.5)

# Here we test all the possible "generic" ways of estimating `TERenyiJizba`.
est_diff = EntropyDecomposition(def, LeonenkoProzantoSavani(Renyi(); k=3))
@test information(est_diff, x, z) isa Real
@test information(est_diff, x, z, y) isa Real

est_disc = EntropyDecomposition(def, PlugIn(Renyi()), ValueBinning(2));
@test information(est_disc, x, z) isa Real
@test information(est_disc, x, z, y) isa Real

# Test `TransferOperator` explicitly
discretization = TransferOperator(RectangularBinning(2))
est_disc = EntropyDecomposition(def, PlugIn(Renyi()), discretization)
@test information(est_disc, x, z) isa Real
@test information(est_disc, x, z, y) isa Real
