using Test
using Associations

def = CMIRenyiSarbu()
@test Associations.min_inputs_vars(def) == 3
@test Associations.max_inputs_vars(def) == 3

# ---------------
# Internals
# ---------------
def = CMIRenyiSarbu()
@test Associations.min_inputs_vars(def) == 3
@test Associations.max_inputs_vars(def) == 3

# Double-sum estimation.
x = rand(["a", "b", "c"], 50)
y = rand(["hello", "yoyo", "heyhey"], 50)
z = rand([1, 2, 5], 50)

# The estimation of probabilities is decoupled from the estimation of the mutual info.
# We could in principle use any probabilities estimator here, but we default to `RelativeAmount`.
p = probabilities(x, y, z) 

@test association(def, p) >= 0.0

est_diff = EntropyDecomposition(def, Kraskov(k=3))
@test_throws ArgumentError association(est_diff, x, z, y)

d = CodifyVariables(OrdinalPatterns(m=3))
est = JointProbabilities(def, d)
@test association(est, x, y, z) isa Real