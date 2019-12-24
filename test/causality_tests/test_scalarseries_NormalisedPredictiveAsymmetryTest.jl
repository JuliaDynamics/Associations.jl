import CausalityToolsBase: RectangularBinning
using StaticArrays 

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

f = 1.0
@test causality(x, y, NormalisedPredictiveAsymmetryTest(test1, f = f)) isa Vector{T} where T
@test causality(x, y, NormalisedPredictiveAsymmetryTest(test2, f = f)) |> typeof <: Real
@test causality(x, y, NormalisedPredictiveAsymmetryTest(test3, f = f)) isa Vector{T} where T
@test causality(x, y, NormalisedPredictiveAsymmetryTest(test4, f = f)) |> typeof<: Real

@test causality(x, y, NormalisedPredictiveAsymmetryTest(test5, f = f)) isa Vector{T} where T
@test causality(x, y, NormalisedPredictiveAsymmetryTest(test6, f = f)) |> typeof<: Real
@test causality(x, y, NormalisedPredictiveAsymmetryTest(test7, f = f)) isa Vector{T} where T
@test causality(x, y, NormalisedPredictiveAsymmetryTest(test8, f = f)) |> typeof<: Real
