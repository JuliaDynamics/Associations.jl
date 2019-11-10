
"""
    EntropyBasedCausalityTest

The supertype of all abstract and composite types representing a causality 
test based on some entropy based measure.

Concrete subtypes are those based on transfer entropy

- [`VisitationFrequencyTest`](@ref)
- [`TransferOperatorGridTest`](@ref)
- [`ExactSimplexIntersectionTest`](@ref)
- [`ApproximateSimplexIntersectionTest`](@ref)

and those that in some manner utilise transfer entropy or other 
information theoretic approaches :

- [`PredictiveAsymmetryTest`](@ref)
"""
abstract type EntropyBasedCausalityTest <: CausalityTest end

export EntropyBasedCausalityTest