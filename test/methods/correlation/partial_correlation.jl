using Test
using StatsBase
x = rand(100)
y = rand(100)
z = rand(100, 2)

@test partial_correlation(x, y, z) ≈ StatsBase.partialcor(x, y, z)
