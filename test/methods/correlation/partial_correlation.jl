using Test
using StatsBase
x = rand(100)
y = rand(100)
z = rand(100, 2)

@test association(PartialCorrelation(), x, y, z) ≈ StatsBase.partialcor(x, y, z)

@test_logs (:warn, "Convenience function `partial_correlation` is deprecated. Use `association(PartialCorrelation(), x, y, z)` instead.") partial_correlation(x, y, z)
