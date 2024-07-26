export TransferEntropyEstimator

"""
The supertype of all dedicated transfer entropy estimators.
"""
abstract type TransferEntropyEstimator{M} <: MultivariateInformationMeasureEstimator{M} end
"""
    estimate_from_marginals(est::TransferEntropyEstimator,
        S::AbstractStateSpaceSet,
        T::AbstractStateSpaceSet,
        T⁺::AbstractStateSpaceSet,
        C::AbstractStateSpaceSet)

Convenience method for [`TransferEntropyEstimator`](@ref)s that allows easier integration
with [`LocalPermutationTest`](@ref). [`Zhu1`](@ref) and [`Lindner`])@ref) uses this method.
"""
function estimate_from_marginals end

# Concrete implementations
include("Zhu1.jl")
include("Lindner.jl")

# convenience
include("Hilbert.jl")
include("SymbolicTransferEntropy.jl")