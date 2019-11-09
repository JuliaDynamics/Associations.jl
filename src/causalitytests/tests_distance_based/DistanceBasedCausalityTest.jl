
"""
    DistanceBasedCausalityTest

The supertype of all abstract and composite types representing a causality 
test based on some sort of distance computation.

Concrete subtypes are 

- [`CrossMappingTest`](@ref)
- [`ConvergentCrossMappingTest`](@ref)
- [`JointDistanceDistributionTest`](@ref)
- [`JointDistanceDistributionTTest`](@ref)
- [`SMeasureTest`](@ref)

"""
abstract type DistanceBasedCausalityTest <: CausalityTest end

export DistanceBasedCausalityTest