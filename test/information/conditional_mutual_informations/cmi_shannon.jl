
using Test
using CausalityTools
using Random
rng = MersenneTwister(1234)
def = CMIShannon()

# ---------------
# Input checks
# ---------------
def = CMIShannon()
@test_throws ArgumentError EntropyDecomposition(def, LeonenkoProzantoSavani(Renyi()))
@test_throws ArgumentError EntropyDecomposition(def, PlugIn(Renyi()), OrdinalPatterns(m=2), RelativeAmount())

# ---------------------------------------------------------------------------------------
# Here we test all the possible "generic" ways of estimating `CMIShannon`.
# Measure-specific tests are in the dedicated estimator test files, e.g. `FPVP.jl`.
# ---------------------------------------------------------------------------------------
# Categorical variables work for `JointProbabilities`
x = rand(["a", "b", "c"], 50)
y = rand(["hello", "yoyo", "heyhey"], 50)
z = rand([1, 2, 5], 50)
est = JointProbabilities(def, UniqueElements())
@test information(est, x, y, z) ≥ 0

x = randn(rng, 50)
y = randn(rng, 50)
z = randn(rng, 50)
est_diff = EntropyDecomposition(def, Kraskov(k=3))
@test information(est_diff, x, z, y) isa Real

est_disc = EntropyDecomposition(def, PlugIn(Shannon()), ValueBinning(2));
@test information(est_disc, x, z, y) isa Real

est_mi = MIDecomposition(def, KSG1())
@test information(est_mi, x, z, y) isa Real

# ---------------
# Pretty printing
# ---------------
out_mi = repr(MIDecomposition(def, KSG1()))
out_hdiff = repr(EntropyDecomposition(def, Kraskov()))
out_hdisc = repr(EntropyDecomposition(def, PlugIn(Shannon()), ValueBinning(2)))

@test occursin("Iₛ(X, Y | Z) = Iₛ(X; Y, Z) + Iₛ(X; Z)", out_mi)
@test occursin("Iₛ(X, Y | Z) = Hₛ(X,Z) + Hₛ(Y,Z) - Hₛ(X,Y,Z) - Hₛ(Z)", out_hdisc)
@test occursin("Iₛ(X, Y | Z) = hₛ(X,Z) + hₛ(Y,Z) - hₛ(X,Y,Z) - hₛ(Z)", out_hdiff)