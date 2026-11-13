using Test
using Associations
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
# JointProbabilities
# ::::::::::::::::::::::::
x = StateSpaceSet(rand(rng, 10, 2)); 
y = StateSpaceSet(rand(rng, 10, 2));
d_row = CodifyPoints(OrdinalPatternEncoding{2}()); 
@test association(JointProbabilities(MIShannon(), d_row), x, y) ≥ 0.0

d_col = CodifyVariables(OrdinalPatterns(m = 2)); 
@test association(JointProbabilities(MIShannon(), d_col), x, y) ≥ 0.0

# ::::::::::::::::::::::::
# Decomposition estimators
# ::::::::::::::::::::::::
x = randn(rng, 50)
y = randn(rng, 50)
est_diff = EntropyDecomposition(def, Kraskov(k=3))
@test association(est_diff, x, y) isa Real

est_disc = EntropyDecomposition(def, PlugIn(Shannon()), CodifyVariables(ValueBinning(2)));
@test association(est_disc, x, y) isa Real

binning = FixedRectangularBinning(0, 1, 3)
disc = CodifyVariables(ValueBinning(binning))
est_bin = EntropyDecomposition(def, PlugIn(Shannon()), disc)
@test association(est_bin, x, y) >= 0.0


# ::::::::::::::::::::::::
# Dedicated estimators
# ::::::::::::::::::::::::
# On vector-valued inputs
def = MIShannon()
x, y = rand(rng, 100), rand(rng, 100)
X, Y = StateSpaceSet(x), StateSpaceSet(y)
@test association(KSG1(def, k = 2), x, y) isa Real
@test association(KSG2(def, k = 2), x, y) isa Real
@test association(GaoOhViswanath(def, k = 2), x, y) isa Real
@test association(GaoKannanOhViswanath(def, k = 2), x, y) isa Real
@test association(GaussianMI(def), x, y) isa Real
@test association(GaussianMI(def), X, Y) isa Real

# On `StateSpaceSet`s
data = [rand(rng, 50, 2) for i = 1:2]
x, y = StateSpaceSet.(data)
def = MIShannon()
@test association(KSG1(def, k = 2), x, y) isa Real
@test association(KSG2(def, k = 2), x, y) isa Real
@test association(GaoOhViswanath(def, k = 2), x, y) isa Real
@test association(GaoKannanOhViswanath(def, k = 2), x, y) isa Real
@test association(GaussianMI(def, normalize = false), x, y) isa Real
@test association(GaussianMI(def; normalize = true), x, y) isa Real

# input is vector + dataset
x = rand(rng, 30)
Y = StateSpaceSet(rand(rng, 30))
@test association(KSG1(def, k = 2), x, Y) isa Real
@test association(KSG2(def, k = 2), x, Y) isa Real
@test association(GaoOhViswanath(def, k = 2), x, Y) isa Real
@test association(GaoKannanOhViswanath(def, k = 2), x, Y) isa Real
@test association(GaussianMI(def, normalize = false), x, Y) isa Real
@test association(GaussianMI(def; normalize = true), x, Y) isa Real

# ---------------
# Pretty printing
# ---------------
def = MIShannon()
out_hdiff = repr(EntropyDecomposition(def, Kraskov()))
out_hdisc = repr(EntropyDecomposition(def, PlugIn(Shannon()), CodifyVariables(ValueBinning(2))))

@test occursin("Iₛ(X, Y) = Hₛ(X) + Hₛ(Y) - Hₛ(X, Y)", out_hdisc)
@test occursin("Iₛ(X, Y) = hₛ(X) + hₛ(Y) - hₛ(X, Y)", out_hdiff)
# ---------------------------------------------------------------
# `GaussianMI` against the exact value, with a multi-column input
# ---------------------------------------------------------------
# a and b are correlated (ρ = 0.8), and y = a + e. Then var(y) = 2 and var(y | a, b) = 1, so
# I((a, b); y) = ½ log 2 nats. (-½ log det of the full correlation matrix would instead give
# the total correlation I((a, b); y) + I(a; b) = ½ log 2 - ½ log(1 - 0.8²) ≈ 0.86 nats.)
let rng = Xoshiro(1234), N = 100_000
    a = randn(rng, N)
    b = 0.8 .* a .+ 0.6 .* randn(rng, N)
    y = a .+ randn(rng, N)
    X = StateSpaceSet(a, b)
    mi_true = 0.5 * log(2)
    for normalize in (true, false)
        est = GaussianMI(MIShannon(; base = ℯ); normalize)
        @test association(est, X, y) ≈ mi_true atol = 0.01
        @test association(est, y, X) ≈ mi_true atol = 0.01
    end
end
