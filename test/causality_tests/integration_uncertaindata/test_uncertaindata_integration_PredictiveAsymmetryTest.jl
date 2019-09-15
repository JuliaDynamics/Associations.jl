
##################################################################################################
# Integration with UncertainData.jl (sampling from full supports of the furnishing distributions)
#################################################################################################

import CausalityToolsBase: RectangularBinning

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
vftest = TransferOperatorGridTest(binning = binnings, ηs = -5:5)
test = PredictiveAsymmetryTest(predictive_test = vftest)

@test causality(x, y, test) isa Array{T, 1} where T
@test causality(uvals_x, uvals_y, test) isa Array{T, 1} where T
@test causality(x, uvals_y, test) isa Array{T, 1} where T
@test causality(uvals_x, y, test) isa Array{T, 1} where T
@test causality(UVX, UVY, test) isa Array{T, 1} where T
@test causality(x, UVY, test) isa Array{T, 1} where T
@test causality(UVX, y, test) isa Array{T, 1} where T



################################################################
# Integration with UncertainData.jl (with sampling constraints)
################################################################
onevar_constraints = ConstrainedResampling(TruncateStd(1))
twovar_constraints = ConstrainedResampling(TruncateStd(2), TruncateStd(1))

test_tog = TransferOperatorGridTest(binning = RectangularBinning(5), ηs = -3:3) 
test_pa = PredictiveAsymmetryTest(predictive_test = test_tog)

@test causality(x, y, twovar_constraints, test_pa) isa Vector{T} where T<:Real
@test causality(uvals_x, uvals_y, twovar_constraints, test_pa) isa Vector{T} where T<:Real
@test causality(x, uvals_y, onevar_constraints, test_pa) isa Vector{T} where T<:Real
@test causality(uvals_x, y, onevar_constraints, test_pa) isa Vector{T} where T<:Real
@test causality(UVX, UVY, twovar_constraints, test_pa) isa Vector{T} where T<:Real
@test causality(uvals_x, UVY, twovar_constraints, test_pa) isa Vector{T} where T<:Real
@test causality(UVX, uvals_y, twovar_constraints, test_pa) isa Vector{T} where T<:Real
@test causality(x, UVY, onevar_constraints, test_pa) isa Vector{T} where T<:Real
@test causality(UVX, y, onevar_constraints, test_pa) isa Vector{T} where T<:Real