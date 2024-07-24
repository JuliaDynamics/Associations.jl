using Test
using CausalityTools

# Double-sum estimation.
x = rand(["a", "b", "c"], 200)
y = rand(["hello", "yoyo", "heyhey"], 200)


# ---------------------------------------------------------------------------------------
# Test all possible ways of estimating `MIRenyiSarbu`.
# ---------------------------------------------------------------------------------------
est = JointProbabilities(MIRenyiSarbu(), UniqueElements())
@test association(est, x, y) ≥ 0.0 # we don't have any better analytical numbers here.

p = probabilities(x, y) 
@test association(MIRenyiSarbu(), p) isa Real # we don't have any better analytical numbers here.
