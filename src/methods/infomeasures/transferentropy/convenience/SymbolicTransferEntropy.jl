export SymbolicTransferEntropy

"""
    SymbolicTransferEntropy <: TransferEntropyEstimator
    SymbolicTransferEntropy(; m = 3, τ = 1, lt = ComplexityMeasures.isless_rand

A convenience estimator for symbolic transfer entropy (Stanieck & Lenertz,
2008)[^Stanieck2008].

## Description

[Symbolic transfer entropy](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.158101)
consists of two simple steps. First, the input time series are embedded with embedding
lag `m` and delay `τ`. The ordinal patterns of the embedding vectors are then encoded
using [`SymbolicPermutation`](@ref) with [`marginal_encodings`](@ref). This transforms the
input time series into integer time series using [`OrdinalPatternEncoding`](@ref).

Transfer entropy is then estimated as usual on the encoded timeseries with
[`transferentropy`](@ref) and the [`CountOccurrences`](@ref) naive frequency estimator.

[^Stanieck2008]:
    Staniek, M., & Lehnertz, K. (2008). Symbolic transfer entropy. Physical review letters,
    100(15), 158101.
"""
Base.@kwdef struct SymbolicTransferEntropy <: TransferEntropyEstimator
    m::Int = 3
    τ::Int = 1
    lt::Function = ComplexityMeasures.isless_rand
end


function transferentropy(measure::TransferEntropy, est::SymbolicTransferEntropy,
    x::AbstractVector...)
    (; m, τ, lt) = est
    est = SymbolicPermutation(; m, τ, lt)
    s = marginal_encodings(est, x...)
    transferentropy(measure, CountOccurrences(), s...)
end
