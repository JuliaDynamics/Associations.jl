import CausalityToolsBase: RectangularBinning

t1 = TransferOperatorGridTest(binning = RectangularBinning(5), ηs = [-4, -3, 3, 5])
t2 = TransferOperatorGridTest(binning = RectangularBinning(5), ηs = 1:5)
@test_throws ArgumentError PredictiveAsymmetryTest(t1)
@test_throws ArgumentError PredictiveAsymmetryTest(t2)

x, y = rand(1000), rand(1000)
single_binning = RectangularBinning(5) 
multiple_binnings = [RectangularBinning(n) for n = 3:5] 

test1 = TransferOperatorGridTest(binning = single_binning, ηs = -10:10)
test2 = TransferOperatorGridTest(binning = single_binning, ηs = 5)
test3 = TransferOperatorGridTest(binning = multiple_binnings, ηs = -10:10)
test4 = TransferOperatorGridTest(binning = multiple_binnings, ηs = 5)

test5 = VisitationFrequencyTest(binning = single_binning, ηs = -10:10)
test6 = VisitationFrequencyTest(binning = single_binning, ηs = 5)
test7 = VisitationFrequencyTest(binning = multiple_binnings, ηs = -10:10)
test8 = VisitationFrequencyTest(binning = multiple_binnings, ηs = 5)


@test causality(x, y, PredictiveAsymmetryTest(test1)) isa Vector{<:Real}
@test causality(x, y, PredictiveAsymmetryTest(test2)) |> typeof <: Real
@test causality(x, y, PredictiveAsymmetryTest(test3)) isa Vector{<:Real}
@test causality(x, y, PredictiveAsymmetryTest(test4)) |> typeof<: Real

@test causality(x, y, PredictiveAsymmetryTest(test5)) isa Vector{<:Real}
@test causality(x, y, PredictiveAsymmetryTest(test6)) |> typeof<: Real
@test causality(x, y, PredictiveAsymmetryTest(test7)) isa Vector{<:Real}
@test causality(x, y, PredictiveAsymmetryTest(test8)) |> typeof<: Real

################################################################
# Integration with UncertainData.jl
################################################################
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