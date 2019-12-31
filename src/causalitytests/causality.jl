

"""
    causality(source::AbstractVector, target::AbstractVector, test::CausalityTest)

Test for a causal influence from `source` to `target` using the provided causality `test`.
    
## Examples

```julia
x, y = rand(300), rand(300)

# Define some causality tests and apply them to `x` and `y`.
test_ccm = ConvergentCrossMappingTest(timeseries_lengths = [45, 50], n_reps = 20)
test_cm = CrossMappingTest(n_reps = 10)
test_vf = VisitationFrequencyTest(binning = RectangularBinning(5), ηs = 1:5)
test_tog = TransferOperatorGridTest(binning = RectangularBinning(5), ηs = 1:5)
test_nn = NearestNeighbourMITest(ηs = 1:5)
test_jdd = JointDistanceDistributionTest()
test_jddt = JointDistanceDistributionTTest()
test_pa = PredictiveAsymmetryTest(
    VisitationFrequencyTest(binning = RectangularBinning(5), ηs = -5:5))
test_pan = NormalisedPredictiveAsymmetryTest(f = 1.0,
    VisitationFrequencyTest(binning = RectangularBinning(5), ηs = -5:5))

causality(x, y, test_ccm)
causality(x, y, test_cm)
causality(x, y, test_vf)
causality(x, y, test_nn)
causality(x, y, test_tog)
causality(x, y, test_jdd)
causality(x, y, test_jddt)
causality(x, y, test_pa)
causality(x, y, test_pan)
```
"""
function causality(source::AbstractVector, target::AbstractVector, test::CausalityTest) end

"""
    causality(source::AbstractVector, target::AbstractVector, cond::AbstractVector, 
        test::CausalityTest)

Test for a causal influence from `source` to `target` conditioned on `cond`, 
using the provided causality `test`. 

Conditional analysis is applicable to the following tests:

- [`VisitationFrequencyTest`](@ref)
- [`TransferOperatorGridTest`](@ref)
- [`PredictiveAsymmetryTest`](@ref)
- [`NormalisedPredictiveAsymmetryTest`](@ref)

## Examples

```julia
x, y, z = rand(300), rand(300), rand(300)

# Transfer entropy from `x` to `y`, conditioned on `z`. Return one 
# transfer entropy value per prediction lag `η`.
test_vf = VisitationFrequencyTest(binning = RectangularBinning(5), ηs = 1:5)

causality(x, y, z, test_vf)
```
"""
function causality(source::AbstractVector, target::AbstractVector, cond::AbstractVector, test::CausalityTest) end


function summarise(test::CausalityTest)
    _type = typeof(test)

    strs = ["$fn = $(test.:($fn))" for fn in fieldnames(_type)]
    return "$_type" * "(" * join(strs, ", ") * ")"
end

Base.show(io::IO, test::CausalityTest) = print(io, summarise(test))
