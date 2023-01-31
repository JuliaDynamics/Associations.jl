using Test
using StatsBase
using Statistics
x = rand(100)
y = rand(100)

pearson_correlation(x, y) ≈ Statistics.cor(x, y)
