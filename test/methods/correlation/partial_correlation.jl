using Test
using StatsBase
using Statistics
x = rand(100)
y = rand(100)
z = rand(100, 2)

partial_correlation(x, y, z) ≈ Statistics.partialcor(x, y, z)
