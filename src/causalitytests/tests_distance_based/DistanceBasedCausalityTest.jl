
"""
    DistanceBasedCausalityTest

The supertype of all abstract and composite types representing a causality 
test based on some sort of distance computation.
"""
abstract type DistanceBasedCausalityTest <: CausalityTest end

export DistanceBasedCausalityTest