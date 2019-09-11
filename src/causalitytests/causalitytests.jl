import CausalityToolsBase: CausalityTest
import UncertainData
import UncertainData: 
    AbstractUncertainValue, 
    AbstractUncertainValueDataset, 
    AbstractUncertainIndexValueDataset,
    resample

"""
    DistanceBasedCausalityTest

The supertype of all abstract and composite types representing a causality 
test based on some sort of distance computation.
"""
abstract type DistanceBasedCausalityTest <: CausalityTest end


"""
    EntropyBasedCausalityTest

The supertype of all abstract and composite types representing a causality 
test based on some entropy based measure.
"""
abstract type EntropyBasedCausalityTest <: CausalityTest end

"""
TransferEntropyTest

The supertype of all abstract and composite types representing a transfer 
entropy causality test.
"""
abstract type TransferEntropyTest <: EntropyBasedCausalityTest end

"""
    directionality(source, target, params::CausalityTest)
    causality(source, target, params::CausalityTest)

Test whether there is a dynamical influence from `source` to `target` by 
applying the causality estimator specified by `params`. 
"""
directionality(source, target, params::CausalityTest)

"""
    causality(source, target, params::CausalityTest)

Test whether there is a dynamical influence from `source` to `target` by 
applying the causality estimator specified by `params`, which must be 
a valid causality test.

## Examples 

```julia
x, y = rand(500), rand(500)

# Test for a causal relation between `x` and `y` using a cross mapping test 
# with default parameters
causality(x, y, CrossMappingTest())
          
# Test for a causal relation between `x` and `y` using a transfer entropy test 
# with default parameters over a rectangular binning with five equidistant 
# intervals along each axis and prediction lag 1
binning =  RectangularBinning(5)
causality(x, y, TransferOperatorGrid(binning = binning, ηs = 1))

# Test for a causal relation between `x` and `y` using a predictive asymmetry 
# test with a transfer entropy test as the underlying test, using 
# default parameters over a rectangular binning with ten equidistant 
# intervals along each axis and prediction lags -10 to 10
binning =  RectangularBinning(10)
te_test = TransferOperatorGridTest(binning = binning, ηs = -10:10)
causality(x, y, PredictiveAsymmetryTest(predictive_test = te_test))
```
"""
causality(source, target, params::CausalityTest)


function summarise(test::CausalityTest)
    _type = typeof(test)
    
    strs = ["$fn = $(test.:($fn))" for fn in fieldnames(_type)]
    return "$_type" * "(" * join(strs, ", ") * ")"
end

Base.show(io::IO, test::CausalityTest) = print(io, summarise(test))

# Resample vectors

UncertainData.resample(v) = v
UncertainData.resample(v::Vector{T}) where T = v

include("distance_based_tests/JointDistanceDistributionTest.jl")
include("distance_based_tests/CrossMappingTest.jl")
include("distance_based_tests/ConvergentCrossMappingTest.jl")
include("entropy_based_tests/VisitationFrequencyTest.jl")
include("entropy_based_tests/TransferOperatorGridTest.jl")
include("predictive_asymmetry/PredictiveAsymmetryTest.jl")


export 
    causality,
    DistanceBasedCausalityTest,
    EntropyBasedCausalityTest,
    TransferEntropyTest,
    resample