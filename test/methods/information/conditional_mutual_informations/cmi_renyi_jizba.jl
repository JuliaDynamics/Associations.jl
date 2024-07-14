using CausalityTools
using Test

# ---------------
# Input checks
# ---------------
def = CMIRenyiJizba()
@test_throws ArgumentError EntropyDecomposition(def, LeonenkoProzantoSavani(Shannon()))
@test_throws ArgumentError EntropyDecomposition(def, PlugIn(Shannon()), OrdinalPatterns(m=2), RelativeAmount())

# ---------------
# Pretty printing
# ---------------
out_hdiff = repr(EntropyDecomposition(def, LeonenkoProzantoSavani(Renyi())))
out_hdisc = repr(EntropyDecomposition(def, PlugIn(Renyi()), ValueBinning(2)))
@test occursin("Iᵣⱼ(X, Y | Z) = Hᵣ(X,Z) + Hᵣ(Y,Z) - Hᵣ(X,Y,Z) - Hᵣ(Z)", out_hdisc)
@test occursin("Iᵣⱼ(X, Y | Z) = hᵣ(X,Z) + hᵣ(Y,Z) - hᵣ(X,Y,Z) - hᵣ(Z)", out_hdiff)


# ---------------------------------------------------------------------------------------
# Test all possible ways of estimating `CMIRenyiJizba`.
# ---------------------------------------------------------------------------------------
x = randn(rng, 50)
y = randn(rng, 50)
z = randn(rng, 50)
est_diff = EntropyDecomposition(def, LeonenkoProzantoSavani(Renyi(), k=3))
@test association(est_diff, x, z, y) isa Real

est_disc = EntropyDecomposition(def, PlugIn(Renyi()), ValueBinning(2));
@test association(est_disc, x, z, y) isa Real

est_joint = JointProbabilities(def, ValueBinning(2))
@test  association(est_joint, x, y, z) isa Real

# ---------------
# Pretty printing
# ---------------
out_hdiff = repr(EntropyDecomposition(def, LeonenkoProzantoSavani(Renyi())))
out_hdisc = repr(EntropyDecomposition(def, PlugIn(Renyi()), ValueBinning(2)))
@test occursin("Iᵣⱼ(X, Y | Z) = Hᵣ(X,Z) + Hᵣ(Y,Z) - Hᵣ(X,Y,Z) - Hᵣ(Z)", out_hdisc)
@test occursin("Iᵣⱼ(X, Y | Z) = hᵣ(X,Z) + hᵣ(Y,Z) - hᵣ(X,Y,Z) - hᵣ(Z)", out_hdiff)