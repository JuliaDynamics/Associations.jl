using Test
using CausalityTools
# ---------------
# Internals
# ---------------
def = CMITsallisPapapetrou()
@test CausalityTools.min_inputs_vars(def) == 3
@test CausalityTools.max_inputs_vars(def) == 3

# Double-sum estimation.
x = rand(["a", "b", "c"], 50)
y = rand(["hello", "yoyo", "heyhey"], 50)
z = rand([1, 2, 5], 50)

# The estimation of probabilities is decoupled from the estimation of the mutual info.
# We could in principle use any probabilities estimator here, but we default to `RelativeAmount`.
p = probabilities(x, y, z) 
@test association(CMITsallisPapapetrou(), p) >= 0.0