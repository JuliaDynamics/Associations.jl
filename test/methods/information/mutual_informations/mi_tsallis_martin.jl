using Test
using CausalityTools
using Random
rng = Xoshiro(1234)

x = rand(rng, ["a", "b", "c"], 200)
y = rand(rng, ["hello", "yoyo", "heyhey"], 200)


# ---------------------------------------------------------------------------------------
# Test all possible ways of estimating `MITsallisMartin`.
# ---------------------------------------------------------------------------------------
def = MITsallisMartin()

# Directly from probabilities
p = probabilities(x, y)
@test association(def, p) isa Real # we don't have any better analytical numbers here.
@test_throws ArgumentError association(MITsallisMartin(q = 1), p)

# `JointProbabilities` estimator
est = JointProbabilities(MITsallisMartin(), UniqueElements())
@test association(est, x, y) isa Real 

# Discrete entropy decomposition
est_disc = EntropyDecomposition(def, PlugIn(Tsallis()), CodifyVariables(UniqueElements()), AddConstant())
@test association(est_disc, x, y) isa Real

# Differential entropy decomposition (on numerical data)
x, y = rand(rng, 100), rand(rng, 100)
est_diff = EntropyDecomposition(def, LeonenkoProzantoSavani(Tsallis(q= 2)))
@test association(est_diff, x, y) isa Real
