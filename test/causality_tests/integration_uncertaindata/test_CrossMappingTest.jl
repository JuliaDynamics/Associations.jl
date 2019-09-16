import Distributions: Normal
import UncertainData: UncertainValue, UncertainValueDataset, resample

# Vectors of uncertain values
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

################################################################
# Integration with UncertainData.jl
################################################################
ctest = CrossMappingTest(n_reps = 10)

# We should get a vector of cross map values
@test causality(x, y, ctest) isa Vector{T} where T
@test causality(uvals_x, uvals_y, ctest) isa Vector{T} where T
@test causality(x, uvals_y, ctest) isa Vector{T} where T
@test causality(uvals_x, y, ctest) isa Vector{T} where T
@test causality(UVX, UVY, ctest) isa Vector{T} where T
@test causality(x, UVY, ctest) isa Vector{T} where T
@test causality(UVX, y, ctest) isa Vector{T} where T

# n_reps = 10, so we should get 10 cross map values
@test causality(x, y, ctest) |> length == 10
@test causality(uvals_x, uvals_y, ctest) |> length == 10
@test causality(x, uvals_y, ctest) |> length == 10
@test causality(uvals_x, y, ctest) |> length == 10
@test causality(UVX, UVY, ctest) |> length == 10
@test causality(x, UVY, ctest) |> length == 10
@test causality(UVX, y, ctest) |> length == 10



################################################################
# Integration with UncertainData.jl (with sampling constraints)
################################################################
onevar_constraints = ConstrainedResampling(TruncateStd(1))
twovar_constraints = ConstrainedResampling(TruncateStd(2), TruncateStd(1))

test_cm = CrossMappingTest() 
@test causality(x, y, test_cm, twovar_constraints) isa Vector{T} where T<:Real
@test causality(uvals_x, uvals_y, test_cm, twovar_constraints) isa Vector{T} where T<:Real
@test causality(x, uvals_y, test_cm, onevar_constraints) isa Vector{T} where T<:Real
@test causality(uvals_x, y, test_cm, onevar_constraints) isa Vector{T} where T<:Real
@test causality(UVX, UVY, test_cm, twovar_constraints) isa Vector{T} where T<:Real
@test causality(uvals_x, UVY, test_cm, twovar_constraints) isa Vector{T} where T<:Real
@test causality(UVX, uvals_y, test_cm, twovar_constraints) isa Vector{T} where T<:Real
@test causality(x, UVY, test_cm, onevar_constraints) isa Vector{T} where T<:Real
@test causality(UVX, y, test_cm, onevar_constraints) isa Vector{T} where T<:Real