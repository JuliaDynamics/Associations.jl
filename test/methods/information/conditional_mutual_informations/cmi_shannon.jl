
using Test
using CausalityTools
using Random
rng = MersenneTwister(1234)

# ---------------
# Internals
# ---------------
def = CMIShannon()
@test CausalityTools.min_inputs_vars(def) == 3
@test CausalityTools.max_inputs_vars(def) == 3

# ---------------
# Input checks
# ---------------
def = CMIShannon()
@test_throws ArgumentError EntropyDecomposition(def, LeonenkoProzantoSavani(Renyi()))
@test_throws ArgumentError EntropyDecomposition(def, PlugIn(Renyi()), OrdinalPatterns(m=2), RelativeAmount())

# ---------------------------------------------------------------------------------------
# Test all possible ways of estimating `CMIShannon`.
# ---------------------------------------------------------------------------------------

# ::::::::::::::::::::::::
# PMF
# ::::::::::::::::::::::::
x = rand(["a", "b", "c"], 50)
y = rand(["hello", "yoyo", "heyhey"], 50)
z = rand([1, 2, 5], 50)
est = JointProbabilities(def, UniqueElements())
@test association(est, x, y, z) ≥ 0

# ::::::::::::::::::::::::
# Decomposition estimators
# ::::::::::::::::::::::::
x = randn(rng, 50)
y = randn(rng, 50)
z = randn(rng, 50)
est_diff = EntropyDecomposition(def, Kraskov(k=3))
@test association(est_diff, x, z, y) isa Real

est_disc = EntropyDecomposition(def, PlugIn(Shannon()), ValueBinning(2));
@test association(est_disc, x, z, y) isa Real

est_mi = MIDecomposition(def, KSG1())
@test association(est_mi, x, z, y) isa Real


# ::::::::::::::::::::::::
# Dedicated estimators
# ::::::::::::::::::::::::
# On vector-valued inputs
def = CMIShannon()
x, y, z = rand(rng, 100), rand(rng, 100), rand(rng, 100);
@test association(MesnerShalizi(def, k = 2), x, y, z) isa Real
@test association(FPVP(def, k = 2), x, y, z) isa Real
@test association(Rahimzamani(def, k = 2), x, y, z) isa Real
@test association(GaussianCMI(def), x, y, z) isa Real

# On `StateSpaceSet`s
data = [rand(rng, 50, 2) for i = 1:3]
x, y, z = StateSpaceSet.(data)
def = CMIShannon()
@test association(MesnerShalizi(def, k = 2), x, y, z) isa Real
@test association(FPVP(def, k = 2), x, y, z) isa Real
@test association(Rahimzamani(def, k = 2), x, y, z) isa Real
@test association(GaussianCMI(def), x, y, z) isa Real


# ---------------
# Pretty printing
# ---------------
out_mi = repr(MIDecomposition(def, KSG1()))
out_hdiff = repr(EntropyDecomposition(def, Kraskov()))
out_hdisc = repr(EntropyDecomposition(def, PlugIn(Shannon()), ValueBinning(2)))

@test occursin("Iₛ(X, Y | Z) = Iₛ(X; Y, Z) + Iₛ(X; Z)", out_mi)
@test occursin("Iₛ(X, Y | Z) = Hₛ(X,Z) + Hₛ(Y,Z) - Hₛ(X,Y,Z) - Hₛ(Z)", out_hdisc)
@test occursin("Iₛ(X, Y | Z) = hₛ(X,Z) + hₛ(Y,Z) - hₛ(X,Y,Z) - hₛ(Z)", out_hdiff)