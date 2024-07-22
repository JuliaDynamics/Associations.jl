using Test
using Statistics
using StateSpaceSets

x = rand(100)
y = rand(100)
X, Y = StateSpaceSet(x), StateSpaceSet(y)

@test association(PearsonCorrelation(), x, y) ≈ Statistics.cor(x, y)
@test association(PearsonCorrelation(), X, Y) ≈ Statistics.cor(x, y)

# Deprecations
@test_logs (:warn, "Convenience function `pearson_correlation` is deprecated. Use `association(PearsonCorrelation(; kwargs...), source, target)` instead.") pearson_correlation(x, y)
