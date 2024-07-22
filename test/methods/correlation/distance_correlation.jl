
using Test
# Analytical tests
# -----------------
a = StateSpaceSet(repeat([1], 100))
# TODO: export this method?
@test CausalityTools.distance_variance(a) == 0.0

v = rand(1000, 3); w = 0.5 .* v .+ 1.2;
@test association(DistanceCorrelation(), v, w) ≈ 1.0
# Comparison with `energy` R package, which is by the authors of the original paper
x = -1.0:0.1:1.0 |> collect
y = map(xᵢ -> xᵢ^3 - 2xᵢ^2 - 3, x)
z = map(yᵢ -> yᵢ^2 - 2yᵢ, y)
dcov = association(DistanceCorrelation(), x, y)
@test round(dcov, digits = 3) == 0.673

# ----------------------------
# Partial distance correlation
# ----------------------------
# Test internals by checking that we get the same answer as in the R `energy` package.
M = reshape([0.0, 0.2, 0.3, 0.2, 0.0, 0.6, 0.3, 0.6, 0.3], 3, 3)
@test CausalityTools.ucenter(M) ≈
    [0 0.15 -0.15;
    0.15 0.0 -0.15;
    -0.15 -0.15 0.0]

@test round(association(DistanceCorrelation(), x, z, y), digits = 5) ≈ round(0.1556139, digits = 5)
@test round(CausalityTools.distance_covariance(x, z, y), digits = 5) ≈ round(0.02379782, digits = 5)

# Deprecations
@test_logs (:warn, "Convenience function `distance_correlation` is deprecated. Use `association(DistanceCorrelation(), x, y)` instead.") distance_correlation(x, y)
@test_logs (:warn, "Convenience function `distance_correlation` is deprecated. Use `association(DistanceCorrelation(), x, y, z)` instead.") distance_correlation(x, y, z)

@test min_inputs_vars(DistanceCorrelation()) == 2
@test max_inputs_vars(DistanceCorrelation()) == 3