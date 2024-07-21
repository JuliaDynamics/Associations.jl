using Test
using CausalityTools

# Double-sum estimation.
x = rand(["a", "b", "c"], 200)
y = rand(["hello", "yoyo", "heyhey"], 200)


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

# ---------------
# Pretty printing
# ---------------
def = MIShannon()
out_hdiff = repr(EntropyDecomposition(def, Kraskov()))
out_hdisc = repr(EntropyDecomposition(def, PlugIn(Shannon()), ValueBinning(2)))

@test occursin("Iₛ(X, Y) = Hₛ(X) + Hₛ(Y) - Hₛ(X, Y)", out_hdisc)
@test occursin("Iₛ(X, Y) = hₛ(X) + hₛ(Y) - hₛ(X, Y)", out_hdiff)