using Test
using CausalityTools

# Double-sum estimation.
x = rand(["a", "b", "c"], 200)
y = rand(["hello", "yoyo", "heyhey"], 200)


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
