using Test
using CausalityTools
using Random
rng = MersenneTwister(1234)
def = MIShannon()

# ---------------------------------------------------------------------------------------
# Test all possible ways of estimating `CMIShannon`.
# ---------------------------------------------------------------------------------------

# ::::::::::::::::::::::::
# PMF
# ::::::::::::::::::::::::
x = rand(rng, ["a", "b", "c"], 200)
y = rand(rng, ["hello", "yoyo", "heyhey"], 200)
est = JointProbabilities(def, UniqueElements())
@test association(est, x, y) ≥ 0

# ::::::::::::::::::::::::
# Decomposition estimators
# ::::::::::::::::::::::::
x = randn(rng, 50)
y = randn(rng, 50)
est_diff = EntropyDecomposition(def, Kraskov(k=3))
@test association(est_diff, x, y) isa Real

est_disc = EntropyDecomposition(def, PlugIn(Shannon()), ValueBinning(2));
@test association(est_disc, x, y) isa Real

# ::::::::::::::::::::::::
# Dedicated estimators
# ::::::::::::::::::::::::
# On vector-valued inputs
def = MIShannon()
x, y = rand(rng, 100), rand(rng, 100)
@test association(KSG1(def, k = 2), x, y) isa Real
@test association(KSG2(def, k = 2), x, y) isa Real
@test association(GaoOhViswanath(def, k = 2), x, y) isa Real
@test association(GaoKannanOhViswanath(def, k = 2), x, y) isa Real
@test association(GaussianMI(def), x, y) isa Real

# On `StateSpaceSet`s
data = [rand(rng, 50, 2) for i = 1:2]
x, y = StateSpaceSet.(data)
def = MIShannon()
@test association(KSG1(def, k = 2), x, y) isa Real
@test association(KSG2(def, k = 2), x, y) isa Real
@test association(GaoOhViswanath(def, k = 2), x, y) isa Real
@test association(GaoKannanOhViswanath(def, k = 2), x, y) isa Real
@test association(GaussianMI(def), x, y) isa Real


# ---------------
# Pretty printing
# ---------------
def = MIShannon()
out_hdiff = repr(EntropyDecomposition(def, Kraskov()))
out_hdisc = repr(EntropyDecomposition(def, PlugIn(Shannon()), ValueBinning(2)))

@test occursin("Iₛ(X, Y) = Hₛ(X) + Hₛ(Y) - Hₛ(X, Y)", out_hdisc)
@test occursin("Iₛ(X, Y) = hₛ(X) + hₛ(Y) - hₛ(X, Y)", out_hdiff)