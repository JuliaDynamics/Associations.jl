using Test
using Associations
using Random
using StableRNGs

rng = StableRNG(123)
sys = system(Logistic4Chain(; rng))
x, y, z, w = columns(first(trajectory(sys, 30, Ttr=10000)))
def = TERenyiJizba(base=3, q=0.5)

# Here we test all the possible "generic" ways of estimating `TERenyiJizba`.
est_diff = EntropyDecomposition(def, LeonenkoProzantoSavani(Renyi(); k=3))
@test association(est_diff, x, z) isa Real
@test association(est_diff, x, z, y) isa Real

est_disc = EntropyDecomposition(def, PlugIn(Renyi()), CodifyVariables(ValueBinning(2)));
@test association(est_disc, x, z) isa Real
@test association(est_disc, x, z, y) isa Real

# Test `TransferOperator` decomposition explicitly, because it has a special implementation
precise = true # precise bin edge
discretization = CodifyVariables(TransferOperator(RectangularBinning(2, precise))) #
est_disc = EntropyDecomposition(TERenyiJizba(), PlugIn(Renyi()), discretization);
@test association(est_disc, x, z) isa Real
@test association(est_disc, x, z, y) isa Real

# Check that in the limit of a lot of points, we roughly get the same answer for transfer 
# operator and regular value binning. 
sys = system(Logistic4Chain(; rng))
x, y, z, w = columns(first(trajectory(sys, 10000, Ttr=10000)))

te_def = TERenyiJizba(base=3, q=0.5)
def_renyi = Renyi()

disc_vf = CodifyVariables(ValueBinning(2))
disc_to = CodifyVariables(TransferOperator(RectangularBinning(2, precise))) #

est_disc_vf = EntropyDecomposition(te_def, PlugIn(def_renyi), disc_vf);
est_disc_to = EntropyDecomposition(te_def, PlugIn(def_renyi), disc_to);
te_vf = association(est_disc_vf, x, z)
te_to = association(est_disc_to, x, z)
# See that values are within 1% of each other
@test in_agreement(te_vf, te_to; agreement_threshold=0.01)

# ---------------
# Pretty printing
# ---------------
te_def = TERenyiJizba(base=3, q=0.5)
out_hdiff = repr(EntropyDecomposition(te_def, LeonenkoProzantoSavani(Renyi())))
out_hdisc = repr(EntropyDecomposition(te_def, PlugIn(Renyi()), CodifyVariables(ValueBinning(2))))

@test occursin("TEᵣⱼ(s → t | c) = hᵣ(t⁺, t⁻,c⁻) - hᵣ(t⁻,c⁻) - hᵣ(t⁺,s⁻,t⁻,c⁻) + hᵣ(s⁻,t⁻,c⁻)", out_hdiff)
@test occursin("TEᵣⱼ(s → t | c) = Hᵣ(t⁺, t⁻,c⁻) - Hᵣ(t⁻,c⁻) - Hᵣ(t⁺,s⁻,t⁻,c⁻) + Hᵣ(s⁻,t⁻,c⁻)", out_hdisc)

