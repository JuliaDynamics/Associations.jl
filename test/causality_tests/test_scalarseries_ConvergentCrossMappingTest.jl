x, y = rand(300), rand(300)

ctest = ConvergentCrossMappingTest(timeseries_lengths = [50, 60, 75, 80], n_reps = 10)
 
# One set of cross mappings per time series length
@test causality(x, y, ctest) isa Vector{Vector{T}} where T

# Four time series lengths, so there should be four vectors of cross map values
@test causality(x, y, ctest) |> length == 4

# There should be ten cross map values per time series length
@test causality(x, y, ctest)[1] |> length == 10