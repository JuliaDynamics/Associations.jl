
"""
    EntropyBasedCausalityTest

The supertype of all abstract and composite types representing a causality 
test based on some entropy based measure.

Concrete subtypes are 

- [`VisitationFrequencyTest`](@ref)
- [`TransferOperatorGridTest`](@ref)
- [`ExactSimplexIntersectionTest`](@ref)
- [`ApproximateSimplexIntersectionTest`](@ref)
- [`PredictiveAsymmetryTest`](@ref)
"""
abstract type EntropyBasedCausalityTest <: CausalityTest end

export EntropyBasedCausalityTest