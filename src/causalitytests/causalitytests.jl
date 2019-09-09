import CausalityToolsBase: CausalityTest


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
    directionality(source, target, params::CausalityTest)
    causality(source, target, params::CausalityTest)

Test whether there is a dynamical influence from `source` to `target` by 
applying the causality estimator specified by `params`. 
"""
causality(source, target, params::CausalityTest)


function summarise(test::CausalityTest)
    _type = typeof(test)
    
    strs = ["$fn = $(test.:($fn))" for fn in fieldnames(_type)]
    return "$_type" * "(" * join(strs, ", ") * ")"
end

Base.show(io::IO, test::CausalityTest) = print(io, summarise(test))


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
    TransferEntropyTest