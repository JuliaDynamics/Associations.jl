using Test
using CausalityTools

# Double-sum estimation.
x = rand(["a", "b", "c"], 200)
y = rand(["hello", "yoyo", "heyhey"], 200)


# ---------------------------------------------------------------------------------------
# Test all possible ways of estimating `MITsallisMartin`.
# ---------------------------------------------------------------------------------------
est = JointProbabilities(MITsallisMartin(), UniqueElements())
@test association(est, x, y) isa Real # we don't have any better analytical numbers here.

# The estimation of probabilities is decoupled from the estimation of the mutual info.
# We could in principle use any probabilities estimator here, but we default to `RelativeAmount`.
p = probabilities(x, y) 
@test association(MITsallisMartin(), p) isa Real # we don't have any better analytical numbers here.
@test_throws ArgumentError association(MITsallisMartin(q = 1), p)