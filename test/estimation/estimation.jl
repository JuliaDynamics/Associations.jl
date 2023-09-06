using Test, CausalityTools
using Random
rng = MersenneTwister(1234)
x, y, z = rand(rng, 50), rand(rng, 50), rand(rng, 50)

# Testing for number of input arguments.
@test_throws ArgumentError estimate(MIShannon(), KSG1(), x)
@test_throws ArgumentError estimate(MIShannon(), KSG1(), x, y, z)
