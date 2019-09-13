import CausalityToolsBase: CausalityTest
import UncertainData
import UncertainData: 
    AbstractUncertainValue, 
    AbstractUncertainValueDataset, 
    AbstractUncertainIndexValueDataset,
    resample


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

################################################################
# Distance based causality tests
################################################################
include("distance_based_tests/DistanceBasedCausalityTest.jl")

# Joint distances tests 
# ---------------------
include("distance_based_tests/JointDistancesCausalityTest.jl")
include("distance_based_tests/JointDistanceDistributionTest.jl")
include("distance_based_tests/JointDistanceDistributionTTest.jl")

# Cross mapping tests 
# ---------------------
include("distance_based_tests/CrossMappingTest.jl")
include("distance_based_tests/ConvergentCrossMappingTest.jl")


################################################################
# Entropy based causality tests
################################################################
include("entropy_based_tests/EntropyBasedCausalityTest.jl")

# Transfer entropy causality tests
# ---------------------------------------
include("entropy_based_tests/TransferEntropyCausalityTest.jl")
include("entropy_based_tests/VisitationFrequencyTest.jl")
include("entropy_based_tests/TransferOperatorGridTest.jl")

################################################################
# Predictive asymmetry causality tests
################################################################
include("predictive_asymmetry/PredictiveAsymmetryTest.jl")


export causality