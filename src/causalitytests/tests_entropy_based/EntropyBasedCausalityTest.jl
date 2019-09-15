
"""
    EntropyBasedCausalityTest

The supertype of all abstract and composite types representing a causality 
test based on some entropy based measure.
"""
abstract type EntropyBasedCausalityTest <: CausalityTest end

export EntropyBasedCausalityTest