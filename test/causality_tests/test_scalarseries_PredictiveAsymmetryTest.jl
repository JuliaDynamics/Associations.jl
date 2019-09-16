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
