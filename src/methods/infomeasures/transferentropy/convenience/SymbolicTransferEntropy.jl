export SymbolicTransferEntropy

"""
    SymbolicTransferEntropy <: TransferEntropyEstimator
    SymbolicTransferEntropy(; m = 3, τ = 1, lt = ComplexityMeasures.isless_rand

A convenience estimator for symbolic transfer entropy [Staniek2008](@cite).

## Description

[Symbolic transfer entropy](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.158101)
consists of two simple steps. First, the input time series are embedded with embedding
lag `m` and delay `τ`. The ordinal patterns of the embedding vectors are then encoded
using [`OrdinalPatterns`](@ref) with [`marginal_encodings`](@ref). This transforms the
input time series into integer time series using [`OrdinalPatternEncoding`](@ref).

Transfer entropy is then estimated as usual on the encoded timeseries with
[`transferentropy`](@ref) and the [`CountOccurrences`](@ref) naive frequency estimator.
"""
Base.@kwdef struct SymbolicTransferEntropy <: TransferEntropyEstimator
    m::Int = 3
    τ::Int = 1
    lt::Function = ComplexityMeasures.isless_rand
end


function estimate(measure::TransferEntropy, est::SymbolicTransferEntropy,
    x::AbstractVector...)
    (; m, τ, lt) = est
    est = OrdinalPatterns(; m, τ, lt)
    s = marginal_encodings(est, x...)
    transferentropy(measure, CountOccurrences(), s...)
end
