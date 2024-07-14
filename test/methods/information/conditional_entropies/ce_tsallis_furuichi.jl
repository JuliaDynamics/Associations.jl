using Test
using CausalityTools

# ---------------
# Internals
# ---------------
def = ConditionalEntropyTsallisFuruichi()
@test CausalityTools.min_inputs_vars(def) == 2
@test CausalityTools.max_inputs_vars(def) == 2

p_nonzeros = Probabilities([0.5 0.5; 0.1 0.1 ])
p_zeros = Probabilities([0.5 0.0; 0.1 0.1])

@test association(ConditionalEntropyTsallisFuruichi(q = 1.5), p_nonzeros) isa Real
@test association(ConditionalEntropyTsallisFuruichi(q = 1.5), p_nonzeros) ≥ 0
@test association(ConditionalEntropyTsallisFuruichi(q = 1.5), p_zeros) |> isnan
@test association(ConditionalEntropyTsallisFuruichi(q = 1.5), p_zeros) |> isnan

@test association(ConditionalEntropyTsallisAbe(q = 1.0), p_zeros) ≥ 0 # shannon
