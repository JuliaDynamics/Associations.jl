using Test
using CausalityTools
using Random


# ---------------
# Internals
# ---------------
def = ConditionalEntropyTsallisAbe()
@test CausalityTools.min_inputs_vars(def) == 2
@test CausalityTools.max_inputs_vars(def) == 2


p_nonzeros = Probabilities([0.5 0.5; 0.1 0.1 ])
p_zeros = Probabilities([0.5 0.0; 0.1 0.1])

@test association(ConditionalEntropyTsallisAbe(q = 1.5), p_nonzeros) isa Real
@test association(ConditionalEntropyTsallisAbe(q = 1.5), p_nonzeros) ≥ 0
@test association(ConditionalEntropyTsallisAbe(q = 1.5), p_zeros) isa Real
@test association(ConditionalEntropyTsallisAbe(q = 1.5), p_zeros) ≥ 0

@test association(ConditionalEntropyTsallisAbe(q = 1.0), p_zeros) ≥ 0 # shannon
