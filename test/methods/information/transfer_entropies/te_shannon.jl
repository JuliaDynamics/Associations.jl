using Test
using CausalityTools
using Random
rng = Xoshiro(1234)

# Double-sum estimation.
x = randn(rng, 50)
y = randn(rng, 50)
z = randn(rng, 50)

# Here we test all the possible "generic" ways of estimating `TEShannon`.
# Remaining tests are in the dedicated estimator test files, e.g. `Zhu1.jl`.
def = TEShannon()
est_diff = EntropyDecomposition(def, Kraskov(k=3))
@test association(est_diff, x, z) isa Real
@test association(est_diff, x, z, y) isa Real

est_disc = EntropyDecomposition(def, PlugIn(Shannon()), CodifyVariables(ValueBinning(2)));
@test association(est_disc, x, z) isa Real
@test association(est_disc, x, z, y) isa Real

est_mi = MIDecomposition(def, KSG1())
@test association(est_mi, x, z) isa Real
@test association(est_mi, x, z, y) isa Real

est_cmi = CMIDecomposition(def, FPVP())
@test association(est_cmi, x, z) isa Real
@test association(est_cmi, x, z, y) isa Real

est_zhu = Zhu1(def, k = 3)
@test association(est_zhu, x, z) isa Real
@test association(est_zhu, x, z, y) isa Real

est_lindner = Lindner(def, k = 3)
@test association(est_lindner, x, z) isa Real
@test association(est_lindner, x, z, y) isa Real


# Test `TransferOperator` decomposition explicitly, because it has a special implementation
precise = true # precise bin edge
discretization = CodifyVariables(TransferOperator(RectangularBinning(2, precise))) #
est_disc = EntropyDecomposition(TEShannon(), PlugIn(Shannon()), discretization);
@test association(est_disc, x, z) isa Real
@test association(est_disc, x, z, y) isa Real

# ---------------
# Pretty printing
# ---------------
out_cmi = repr(CMIDecomposition(def, FPVP()))
out_mi = repr(MIDecomposition(def, KSG1()))
out_hdiff = repr(EntropyDecomposition(def, Kraskov()))
out_hdisc = repr(EntropyDecomposition(def, PlugIn(Shannon()), CodifyVariables(ValueBinning(2))))

@test occursin("TEₛ(s → t | c) = Iₛ(t⁺; s⁻ | t⁻, c⁻)", out_cmi)
@test occursin("TEₛ(s → t | c) = Iₛ(t⁺; s⁻, t⁻, c⁻) - Iₛ(t⁺; t⁻, c⁻)", out_mi)
@test occursin("TEₛ(s → t | c) = hₛ(t⁺, t⁻,c⁻) - hₛ(t⁻,c⁻) - hₛ(t⁺,s⁻,t⁻,c⁻) + hₛ(s⁻,t⁻,c⁻)", out_hdiff)
@test occursin("TEₛ(s → t | c) = Hₛ(t⁺, t⁻,c⁻) - Hₛ(t⁻,c⁻) - Hₛ(t⁺,s⁻,t⁻,c⁻) + Hₛ(s⁻,t⁻,c⁻)", out_hdisc)