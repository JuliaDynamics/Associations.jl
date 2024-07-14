using Test
using CausalityTools
using Random
rng = MersenneTwister(1234)

def = MIRenyiJizba(q = 0.5)

# ---------------------------------------------------------------------------------------
# Test all possible ways of estimating `MIRenyiJizba`.
# ---------------------------------------------------------------------------------------
# ::::::::::::::::::::::::
# PMF
# ::::::::::::::::::::::::
x = rand(rng, ["a", "b", "c"], 200);
y = rand(rng, ["hello", "yoyo", "heyhey"], 200);
est = JointProbabilities(def, UniqueElements())
@test association(est, x, y) ≥ 0

# ::::::::::::::::::::::::
# Decomposition estimators
# ::::::::::::::::::::::::
x = randn(rng, 50);
y = randn(rng, 50);
est_diff = EntropyDecomposition(def, LeonenkoProzantoSavani(Renyi(), k=3))
@test association(est_diff, x, y) isa Real

est_disc = EntropyDecomposition(def, PlugIn(Renyi()), ValueBinning(2));
@test association(est_disc, x, y) isa Real

# ---------------
# Pretty printing
# ---------------
def = MIRenyiJizba()
out_hdiff = repr(EntropyDecomposition(def, LeonenkoProzantoSavani(Renyi())))
out_hdisc = repr(EntropyDecomposition(def, PlugIn(Renyi()), ValueBinning(2)))
@test occursin("Iᵣⱼ(X, Y) = Hᵣ(X) + Hᵣ(Y) - Hᵣ(X, Y)", out_hdisc)
@test occursin("Iᵣⱼ(X, Y) = hᵣ(X) + hᵣ(Y) - hᵣ(X, Y)", out_hdiff)