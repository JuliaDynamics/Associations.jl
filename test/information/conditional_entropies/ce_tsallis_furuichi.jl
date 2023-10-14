using Test
using CausalityTools

p_nonzeros = Probabilities([0.5 0.5; 0.1 0.1 ])
p_zeros = Probabilities([0.5 0.0; 0.1 0.1])

@test information(ConditionalEntropyTsallisFuruichi(q = 1.5), p_nonzeros) isa Real
@test information(ConditionalEntropyTsallisFuruichi(q = 1.5), p_nonzeros) ≥ 0
@test information(ConditionalEntropyTsallisFuruichi(q = 1.5), p_zeros) |> isnan
@test information(ConditionalEntropyTsallisFuruichi(q = 1.5), p_zeros) |> isnan
