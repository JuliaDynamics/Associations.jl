using Test
using Statistics
x = rand(100)
y = rand(100)

@test association(PearsonCorrelation(), x, y) ≈ Statistics.cor(x, y)

# Deprecations
@test_logs (:warn, "Convenience function `pearson_correlation` is deprecated. Use `association(PearsonCorrelation(; kwargs...), source, target)` instead.") pearson_correlation(x, y)
