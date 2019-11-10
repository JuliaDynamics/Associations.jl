
"""
    TransferEntropyCausalityTest

The supertype of all abstract and composite types representing a transfer 
entropy causality test.

Concrete subtypes are 

- [`TransferOperatorGridTest`](@ref)
- [`VisitationFrequencyTest`](@ref)
- [`ApproximateSimplexIntersectionTest`](@ref)
- [`ExactSimplexIntersectionTest`](@ref)
"""
abstract type TransferEntropyCausalityTest <: EntropyBasedCausalityTest end

export TransferEntropyCausalityTest