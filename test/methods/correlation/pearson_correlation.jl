using Test
using Statistics
x = rand(100)
y = rand(100)

@test pearson_correlation(x, y) ≈ Statistics.cor(x, y)
