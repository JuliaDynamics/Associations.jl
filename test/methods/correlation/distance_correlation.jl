
using Test
# Analytical tests
# -----------------
a = Dataset(repeat([1], 100))
# TODO: export this method?
@test CausalityTools.distance_variance(a) == 0.0

v = rand(1000, 3); w = 0.5 .* v .+ 1.2;
@test distance_correlation(v, w) ≈ 1.0
# Comparison with `energy` R package, which is by the authors of the original paper
x = -1.0:0.1:1.0 |> collect
y = map(xᵢ -> xᵢ^3 - 2xᵢ^2 - 3, x)
dcov = distance_correlation(x, y)
@test round(dcov, digits = 3) == 0.673
