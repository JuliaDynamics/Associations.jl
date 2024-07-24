export TransferEntropyEstimator

"""
The supertype of all dedicated transfer entropy estimators.
"""
abstract type TransferEntropyEstimator{M} <: MultivariateInformationMeasureEstimator{M} end

# Concrete implementations
include("Zhu1.jl")
include("Lindner.jl")

# convenience
include("Hilbert.jl")
include("SymbolicTransferEntropy.jl")