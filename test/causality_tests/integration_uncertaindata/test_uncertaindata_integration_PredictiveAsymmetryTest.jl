
##################################################################################################
# Integration with UncertainData.jl (sampling from full supports of the furnishing distributions)
#################################################################################################
using StaticArrays
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

@test causality(x, y, test) isa SVector{5, T} where T
@test causality(uvals_x, uvals_y, test) isa SVector{5, T} where T
@test causality(x, uvals_y, test) isa SVector{5, T} where T
@test causality(uvals_x, y, test) isa SVector{5, T} where T
@test causality(UVX, UVY, test) isa SVector{5, T} where T
@test causality(x, UVY, test) isa SVector{5, T} where T
@test causality(UVX, y, test) isa SVector{5, T} where T



################################################################
# Integration with UncertainData.jl (with sampling constraints)
################################################################
onevar_constraints = ConstrainedResampling(TruncateStd(1))
twovar_constraints = ConstrainedResampling(TruncateStd(2), TruncateStd(1))

test_tog = TransferOperatorGridTest(binning = RectangularBinning(5), ηs = -3:3) 
test_pa = PredictiveAsymmetryTest(predictive_test = test_tog)

@test causality(x, y, test_pa, twovar_constraints) isa SVector{3, T} where T
@test causality(uvals_x, uvals_y, test_pa, twovar_constraints) isa SVector{3, T} where T
@test causality(x, uvals_y, test_pa, onevar_constraints) isa SVector{3, T} where T
@test causality(uvals_x, y, test_pa, onevar_constraints) isa SVector{3, T} where T
@test causality(UVX, UVY, test_pa, twovar_constraints) isa SVector{3, T} where T
@test causality(uvals_x, UVY, test_pa, twovar_constraints) isa SVector{3, T} where T
@test causality(UVX, uvals_y, test_pa, twovar_constraints) isa SVector{3, T} where T
@test causality(x, UVY, test_pa, onevar_constraints) isa SVector{3, T} where T
@test causality(UVX, y, test_pa, onevar_constraints) isa SVector{3, T} where T