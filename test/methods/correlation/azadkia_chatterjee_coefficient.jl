using Test
using Random
rng = Xoshiro(12345)
n = 15
x = rand(rng, n)
y = rand(rng, n)
z = rand(rng, n)
m = AzadkiaChatterjeeCoefficient()
@test association(m, x, y) isa Real
@test association(m, x, y, z) isa Real