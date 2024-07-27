using Test
using Associations
using Random
rng = Xoshiro(1234)

x = rand(rng, ["a", "b", "c"], 200)
y = rand(rng, ["hello", "yoyo", "heyhey"], 200)

# ---------------------------------------------------------------------------------------
# Test all possible ways of estimating `MITsallisFuruichi`.
# ---------------------------------------------------------------------------------------
def = MITsallisFuruichi()
# Directly from probabilities
p = probabilities(x, y)
@test association(def, p) ≥ 0.0

# `JointProbabilities` estimator
est = JointProbabilities(def, UniqueElements())
@test association(est, x, y) ≥ 0.0 # we don't have any better analytical numbers here.

# Discrete entropy decomposition (on numerical data)
x, y = rand(rng, 100), rand(rng, 100)
est_disc = EntropyDecomposition(def, PlugIn(Tsallis()), CodifyVariables(OrdinalPatterns()), AddConstant())
@test association(est_disc, x, y) isa Real

# Differential entropy decomposition (on numerical data)
x, y = rand(rng, 100), rand(rng, 100)
est_diff = EntropyDecomposition(def, LeonenkoProzantoSavani(Tsallis(q= 2)))
@test association(est_diff, x, y) isa Real

# ---------------
# Pretty printing
# ---------------
def = MIShannon()
out_hdiff = repr(EntropyDecomposition(def, Kraskov()))
out_hdisc = repr(EntropyDecomposition(def, PlugIn(Shannon()), CodifyVariables(ValueBinning(2))))

@test occursin("Iₛ(X, Y) = Hₛ(X) + Hₛ(Y) - Hₛ(X, Y)", out_hdisc)
@test occursin("Iₛ(X, Y) = hₛ(X) + hₛ(Y) - hₛ(X, Y)", out_hdiff)