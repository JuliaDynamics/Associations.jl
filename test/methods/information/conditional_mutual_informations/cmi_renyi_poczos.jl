
using Test
using Associations
using Random
rng = MersenneTwister(1234)

x = randn(rng, 50)
y = randn(rng, 50)
z = randn(rng, 50)


# ---------------
# Internals
# ---------------
def = CMIRenyiPoczos()
@test Associations.min_inputs_vars(def) == 3
@test Associations.max_inputs_vars(def) == 3

# ---------------------------------------------------------------------------------------
# Test all possible ways of estimating `CMIRenyiPoczos`.
# ---------------------------------------------------------------------------------------
def = CMIRenyiPoczos()
@test association(PoczosSchneiderCMI(def, k = 2), x, y, z) isa Real

data = [rand(rng, 50, 2) for i = 1:3]
x, y, z = StateSpaceSet.(data)
@test association(PoczosSchneiderCMI(def, k = 2), x, y, z) isa Real
