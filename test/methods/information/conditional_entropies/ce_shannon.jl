using Test
using CausalityTools
using Random


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
