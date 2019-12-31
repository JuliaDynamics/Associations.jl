
"""
    TransferEntropyCausalityTest{N}

The supertype of all abstract and composite types representing a transfer 
entropy causality test applied to `N` different prediction lags.

Concrete subtypes are 

- [`TransferOperatorGridTest`](@ref)
- [`VisitationFrequencyTest`](@ref)
- [`ApproximateSimplexIntersectionTest`](@ref)
- [`ExactSimplexIntersectionTest`](@ref)
"""
abstract type TransferEntropyCausalityTest{N} <: EntropyBasedCausalityTest end

export TransferEntropyCausalityTest