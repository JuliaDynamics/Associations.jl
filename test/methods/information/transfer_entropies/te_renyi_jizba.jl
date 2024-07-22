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
@test association(est_diff, x, z) isa Real
@test association(est_diff, x, z, y) isa Real

est_disc = EntropyDecomposition(def, PlugIn(Renyi()), CodifyVariables(ValueBinning(2)));
@test association(est_disc, x, z) isa Real
@test association(est_disc, x, z, y) isa Real

# TODO: how to incorporate TransferOperator estimation explicitly? 
# We probabily need a dedicated estimator.
# # Test `TransferOperator` explicitly
# discretization = CodifyVariables(TransferOperator(RectangularBinning(2, true)))
# est_disc = EntropyDecomposition(def, PlugIn(Renyi()), discretization)
# @test association(est_disc, x, z) isa Real
# @test association(est_disc, x, z, y) isa Real



# ---------------
# Pretty printing
# ---------------
out_hdiff = repr(EntropyDecomposition(def, LeonenkoProzantoSavani(Renyi())))
out_hdisc = repr(EntropyDecomposition(def, PlugIn(Renyi()), CodifyVariables(ValueBinning(2))))

@test occursin("TEᵣⱼ(s → t | c) = hᵣ(t⁺, t⁻,c⁻) - hᵣ(t⁻,c⁻) - hᵣ(t⁺,s⁻,t⁻,c⁻) + hᵣ(s⁻,t⁻,c⁻)", out_hdiff)
@test occursin("TEᵣⱼ(s → t | c) = Hᵣ(t⁺, t⁻,c⁻) - Hᵣ(t⁻,c⁻) - Hᵣ(t⁺,s⁻,t⁻,c⁻) + Hᵣ(s⁻,t⁻,c⁻)", out_hdisc)