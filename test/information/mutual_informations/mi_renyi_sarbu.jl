using Test
using CausalityTools

# Double-sum estimation.
x = rand(["a", "b", "c"], 200)
y = rand(["hello", "yoyo", "heyhey"], 200)

# The estimation of probabilities is decoupled from the estimation of the mutual info.
# We could in principle use any probabilities estimator here, but we default to `RelativeAmount`.
p = probabilities(x, y) 
@test information(MIRenyiSarbu(), p) isa Real # we don't have any better analytical numbers here.
