using Test
using CausalityTools
using Random

p_nonzeros = Probabilities([0.5 0.5; 0.1 0.1 ])
p_zeros = Probabilities([0.5 0.0; 0.1 0.1])

@test information(CEShannon(), p_nonzeros) isa Real
@test information(CEShannon(), p_nonzeros) ≥ 0
@test information(CEShannon(), p_zeros) isa Real
@test information(CEShannon(), p_zeros) ≥ 0
