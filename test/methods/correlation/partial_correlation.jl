using Test
using StatsBase: partialcor
x = rand(100)
y = rand(100)
z = rand(100, 2)

@test association(PartialCorrelation(), x, y, z) ≈ partialcor(x, y, z)

@test_logs (:warn, "Convenience function `partial_correlation` is deprecated. Use `association(PartialCorrelation(), x, y, z)` instead.") partial_correlation(x, y, z)

@test CausalityTools.min_inputs_vars(PartialCorrelation()) == 3
@test CausalityTools.max_inputs_vars(PartialCorrelation()) == Inf