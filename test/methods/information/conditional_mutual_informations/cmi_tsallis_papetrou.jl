using Test
using Associations
# ---------------
# Internals
# ---------------
def = CMITsallisPapapetrou()
@test Associations.min_inputs_vars(def) == 3
@test Associations.max_inputs_vars(def) == 3

# ---------------------------------------------------------------------------------------
# Test all possible ways of estimating `CMIRenyiJizba`.
# ---------------------------------------------------------------------------------------

x = rand(["a", "b", "c"], 50)
y = rand(["hello", "yoyo", "heyhey"], 50)
z = rand([1, 2, 5], 50)

# From raw probabilities
p = probabilities(x, y, z) 
@test association(CMITsallisPapapetrou(), p) >= 0.0

# `JointProbabilities` estimator
d = CodifyVariables(UniqueElements())
est_joint = JointProbabilities(def, d)
@test  association(est_joint, x, y, z) isa Real
