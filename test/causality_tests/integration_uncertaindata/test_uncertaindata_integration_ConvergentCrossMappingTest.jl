
##################################################################################################
# Integration with UncertainData.jl (sampling from full supports of the furnishing distributions)
#################################################################################################

#Vectors of uncertain values
uvals_x = [UncertainValue(Normal, rand(Normal(0, 5)), abs(rand(Normal(0, 3)))) for i = 1:100]
uvals_y = [UncertainValue(Normal, rand(Normal(0, 5)), abs(rand(Normal(0, 3)))) for i = 1:100];

# UncertainValueDataset
UVX = UncertainValueDataset(uvals_x)
UVY = UncertainValueDataset(uvals_y)

# UncertainIndexDataset
UVX_idx = UncertainIndexDataset(uvals_x)
UVY_idx = UncertainIndexDataset(uvals_y)

# Real-valued vectors
x = resample.(uvals_x);
y = resample.(uvals_y);

ctest = ConvergentCrossMappingTest(timeseries_lengths = [50, 60, 75, 80], n_reps = 10)
 
# One set of cross mappings per time series length
@test causality(x, y, ctest) isa Vector{Vector{T}} where T
@test causality(uvals_x, uvals_y, ctest) isa Vector{Vector{T}} where T
@test causality(x, uvals_y, ctest) isa Vector{Vector{T}} where T
@test causality(uvals_x, y, ctest) isa Vector{Vector{T}} where T
@test causality(UVX, UVY, ctest) isa Vector{Vector{T}} where T
@test causality(x, UVY, ctest) isa Vector{Vector{T}} where T
@test causality(UVX, y, ctest) isa Vector{Vector{T}} where T

# Four time series lengths, so there should be four vectors of cross map values
@test causality(x, y, ctest) |> length == 4
@test causality(uvals_x, uvals_y, ctest) |> length == 4
@test causality(x, uvals_y, ctest) |> length == 4
@test causality(uvals_x, y, ctest) |> length == 4
@test causality(UVX, UVY, ctest) |> length == 4
@test causality(x, UVY, ctest) |> length == 4
@test causality(UVX, y, ctest) |> length == 4

# There should be ten cross map values per time series length
@test causality(x, y, ctest)[1] |> length == 10
@test causality(uvals_x, uvals_y, ctest)[1] |> length == 10
@test causality(x, uvals_y, ctest)[1] |> length == 10
@test causality(uvals_x, y, ctest)[1] |> length == 10
@test causality(UVX, UVY, ctest)[1] |> length == 10
@test causality(x, UVY, ctest)[1] |> length == 10
@test causality(UVX, y, ctest)[1] |> length == 10


################################################################
# Integration with UncertainData.jl (with sampling constraints)
################################################################
onevar_constraints = ConstrainedResampling(TruncateStd(1))
twovar_constraints = ConstrainedResampling(TruncateStd(2), TruncateStd(1))

test_ccm = ConvergentCrossMappingTest(timeseries_lengths = [50, 60], n_reps = 5)
@test causality(x, y, twovar_constraints, test_ccm) isa Vector{Vector{T}} where T<:Real
@test causality(uvals_x, uvals_y, twovar_constraints, test_ccm) isa Vector{Vector{T}} where T<:Real
@test causality(x, uvals_y, onevar_constraints, test_ccm) isa Vector{Vector{T}} where T<:Real
@test causality(uvals_x, y, onevar_constraints, test_ccm) isa Vector{Vector{T}} where T<:Real
@test causality(UVX, UVY, twovar_constraints, test_ccm) isa Vector{Vector{T}} where T<:Real
@test causality(uvals_x, UVY, twovar_constraints, test_ccm) isa Vector{Vector{T}} where T<:Real
@test causality(UVX, uvals_y, twovar_constraints, test_ccm) isa Vector{Vector{T}} where T<:Real
@test causality(x, UVY, onevar_constraints, test_ccm) isa Vector{Vector{T}} where T<:Real
@test causality(UVX, y, onevar_constraints, test_ccm) isa Vector{Vector{T}} where T<:Real