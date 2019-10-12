
##################################################################################################
# Integration with UncertainData.jl (sampling from full supports of the furnishing distributions)
#################################################################################################

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

binnings = [RectangularBinning(i) for i = 3:6]
vftest = TransferOperatorGridTest(binning = binnings, ηs = 1)

@test causality(x, y, vftest) isa Array{<:Real, 0}
@test causality(uvals_x, uvals_y, vftest) isa Array{<:Real, 0}
@test causality(x, uvals_y, vftest) isa Array{<:Real, 0}
@test causality(uvals_x, y, vftest) isa Array{<:Real, 0}
@test causality(UVX, UVY, vftest) isa Array{<:Real, 0}
@test causality(x, UVY, vftest) isa Array{<:Real, 0}
@test causality(UVX, y, vftest) isa Array{<:Real, 0}

vftest = TransferOperatorGridTest(binning = binnings, ηs = 1:5)

@test causality(x, y, vftest) isa Array{T, 1} where T
@test causality(uvals_x, uvals_y, vftest) isa Array{T, 1} where T
@test causality(x, uvals_y, vftest) isa Array{T, 1} where T
@test causality(uvals_x, y, vftest) isa Array{T, 1} where T
@test causality(UVX, UVY, vftest) isa Array{T, 1} where T
@test causality(x, UVY, vftest) isa Array{T, 1} where T
@test causality(UVX, y, vftest) isa Array{T, 1} where T

@test causality(x, y, vftest) |> length == 5
@test causality(uvals_x, uvals_y, vftest) |> length == 5
@test causality(x, uvals_y, vftest) |> length == 5
@test causality(uvals_x, y, vftest) |> length == 5
@test causality(UVX, UVY, vftest) |> length == 5
@test causality(x, UVY, vftest) |> length == 5
@test causality(UVX, y, vftest) |> length == 5



################################################################
# Integration with UncertainData.jl (with sampling constraints)
################################################################
onevar_constraints = ConstrainedResampling(TruncateStd(1))
twovar_constraints = ConstrainedResampling(TruncateStd(2), TruncateStd(1))

test_tog = TransferOperatorGridTest(binning = RectangularBinning(5), ηs = 1) 
@test causality(x, y, test_tog, twovar_constraints) isa Array{<:Real, 0}
@test causality(uvals_x, uvals_y, test_tog, twovar_constraints) isa Array{<:Real, 0}
@test causality(x, uvals_y, test_tog, onevar_constraints) isa Array{<:Real, 0}
@test causality(uvals_x, y, test_tog, onevar_constraints) isa Array{<:Real, 0}
@test causality(UVX, UVY, test_tog, twovar_constraints) isa Array{<:Real, 0}
@test causality(uvals_x, UVY, test_tog, twovar_constraints) isa Array{<:Real, 0}
@test causality(UVX, uvals_y, test_tog, twovar_constraints) isa Array{<:Real, 0}
@test causality(x, UVY, test_tog, onevar_constraints) isa Array{<:Real, 0}
@test causality(UVX, y, test_tog, onevar_constraints) isa Array{<:Real, 0}

test_tog = TransferOperatorGridTest(binning = RectangularBinning(5), ηs = -3:3) 
@test causality(x, y, test_tog, twovar_constraints) isa Array{<:Real, 1}
@test causality(uvals_x, uvals_y, test_tog, twovar_constraints) isa Array{<:Real, 1}
@test causality(x, uvals_y, test_tog, onevar_constraints) isa Array{<:Real, 1}
@test causality(uvals_x, y, test_tog, onevar_constraints) isa Array{<:Real, 1}
@test causality(UVX, UVY, test_tog, twovar_constraints) isa Array{<:Real, 1}
@test causality(uvals_x, UVY, test_tog, twovar_constraints) isa Array{<:Real, 1}
@test causality(UVX, uvals_y, test_tog, twovar_constraints) isa Array{<:Real, 1}
@test causality(x, UVY, test_tog, onevar_constraints) isa Array{<:Real, 1}
@test causality(UVX, y, test_tog, onevar_constraints) isa Array{<:Real, 1}