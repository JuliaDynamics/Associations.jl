
##################################################################################################
# Integration with UncertainData.jl (sampling from full supports of the furnishing distributions)
#################################################################################################

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

jtest = JointDistanceDistributionTest(B = 20)

@test causality(x, y, jtest) isa Vector{T} where T
@test causality(uvals_x, uvals_y, jtest) isa Vector{T} where T
@test causality(x, uvals_y, jtest) isa Vector{T} where T
@test causality(uvals_x, y, jtest) isa Vector{T} where T
@test causality(UVX, UVY, jtest) isa Vector{T} where T
@test causality(x, UVY, jtest) isa Vector{T} where T
@test causality(UVX, y, jtest) isa Vector{T} where T

# 20 intervals, so 20 values should be in the joint distance distribution
@test causality(x, y, jtest) |> length == 20
@test causality(uvals_x, uvals_y, jtest) |> length == 20
@test causality(x, uvals_y, jtest) |> length == 20
@test causality(uvals_x, y, jtest) |> length == 20
@test causality(UVX, UVY, jtest) |> length == 20
@test causality(x, UVY, jtest) |> length == 20
@test causality(UVX, y, jtest) |> length == 20

################################################################
# Integration with UncertainData.jl (with sampling constraints)
################################################################
onevar_constraints = ConstrainedResampling(TruncateStd(1))
twovar_constraints = ConstrainedResampling(TruncateStd(2), TruncateStd(1))

test_jdd = JointDistanceDistributionTest() 
@test causality(x, y, twovar_constraints, test_jdd) isa Vector{T} where T<:Real
@test causality(uvals_x, uvals_y, twovar_constraints, test_jdd) isa Vector{T} where T<:Real
@test causality(x, uvals_y, onevar_constraints, test_jdd) isa Vector{T} where T<:Real
@test causality(uvals_x, y, onevar_constraints, test_jdd) isa Vector{T} where T<:Real
@test causality(UVX, UVY, twovar_constraints, test_jdd) isa Vector{T} where T<:Real
@test causality(uvals_x, UVY, twovar_constraints, test_jdd) isa Vector{T} where T<:Real
@test causality(UVX, uvals_y, twovar_constraints, test_jdd) isa Vector{T} where T<:Real
@test causality(x, UVY, onevar_constraints, test_jdd) isa Vector{T} where T<:Real
@test causality(UVX, y, onevar_constraints, test_jdd) isa Vector{T} where T<:Real