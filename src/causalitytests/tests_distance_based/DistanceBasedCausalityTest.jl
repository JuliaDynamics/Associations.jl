
"""
    DistanceBasedCausalityTest{N}

The supertype of all abstract and composite types representing a causality 
test based on some sort of distance computation. 

The type parameter `N` indicates the number of returned elements when 
applying the test.

Concrete subtypes are 

- [`CrossMappingTest`](@ref)
- [`ConvergentCrossMappingTest`](@ref)
- [`JointDistanceDistributionTest`](@ref)
- [`JointDistanceDistributionTTest`](@ref)
- [`SMeasureTest`](@ref)

"""
abstract type DistanceBasedCausalityTest{N} <: CausalityTest end

export DistanceBasedCausalityTest