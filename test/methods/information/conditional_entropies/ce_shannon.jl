using Test
using CausalityTools
using Random
rng = Xoshiro(1234)

# ---------------
# Internals
# ---------------
def = ConditionalEntropyShannon()
@test CausalityTools.min_inputs_vars(def) == 2
@test CausalityTools.max_inputs_vars(def) == 2

p_nonzeros = Probabilities([0.5 0.5; 0.1 0.1 ])
p_zeros = Probabilities([0.5 0.0; 0.1 0.1])

@test association(ConditionalEntropyShannon(), p_nonzeros) isa Real
@test association(ConditionalEntropyShannon(), p_nonzeros) ≥ 0
@test association(ConditionalEntropyShannon(), p_zeros) isa Real
@test association(ConditionalEntropyShannon(), p_zeros) ≥ 0


# ---------------------------------------------------------------------------------------
# Test all possible ways of estimating `ConditionalEntropyShannon`.
# ---------------------------------------------------------------------------------------
# `JointProbabilities` with ` CodifyPoints`
x, y, z = rand(rng, 1:5, 100), rand(rng, 1:5, 100), rand(rng, 1:3, 100)
X = StateSpaceSet(x, z)
Y = StateSpaceSet(y, z)
disc = CodifyPoints(UniqueElementsEncoding(X), UniqueElementsEncoding(Y));
est = JointProbabilities(ConditionalEntropyShannon(), disc);
association(est, X, Y)