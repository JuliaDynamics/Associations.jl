x, y = rand(300), rand(300)

ctest = CrossMappingTest(n_reps = 10)

# We should get a vector of cross map values
@test causality(x, y, ctest) isa Vector{T} where T