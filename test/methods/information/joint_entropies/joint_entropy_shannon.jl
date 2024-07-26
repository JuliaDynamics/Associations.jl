using Test
using CausalityTools

x = rand(["a", "b", "c"], 200)
y = rand(["hello", "yoyo", "heyhey"], 200)

# The estimation of probabilities is decoupled from the estimation of the mutual info.
# We could in principle use any probabilities estimator here, but we default to `RelativeAmount`.
p = probabilities(x, y) 
@test association(JointEntropyShannon(), p) ≥ 0

# TODO: estimation using discretizations..