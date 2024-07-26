using Test
using CausalityTools

# Double-sum estimation.
x = rand(["a", "b", "c"], 200)
y = rand(["hello", "yoyo", "heyhey"], 200)


# ---------------------------------------------------------------------------------------
# Test all possible ways of estimating `MIRenyiSarbu`.
# ---------------------------------------------------------------------------------------

# ::::::::::::::::::::::::
# From raw probabilities
# ::::::::::::::::::::::::
p = probabilities(x, y) 
@test association(MIRenyiSarbu(), p) isa Real # we don't have any better analytical numbers here.

# ::::::::::::::::::::::::
# JointProbabilities
# ::::::::::::::::::::::::
est = JointProbabilities(MIRenyiSarbu(), CodifyVariables(UniqueElements()))
@test association(est, x, y) ≥ 0.0 # we don't have any better analytical numbers here.

x = StateSpaceSet(rand(rng, 10, 2)); 
y = StateSpaceSet(rand(rng, 10, 2));
d_row = CodifyPoints(OrdinalPatternEncoding{2}()); 
@test association(JointProbabilities(MIShannon(), d_row), x, y) ≥ 0.0

