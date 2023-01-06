import ComplexityMeasures: symbolize_for_dispersion
export marginal_encodings

"""
    marginal_encodings(est::ProbabilitiesEstimator, x::AbstractVector...)

Encode/discretize each input vector `xᵢ ∈ x` according to a procedure determined by `est`.

This is useful for computing any discrete information theoretic quantity, and is
used internally by [`contingency_matrix`](@ref).

## Supported estimators

- [`ValueHistogram`](@ref) with [`FixedRectangularBinning`](@ref). Bin visitation
    frequencies are counted in the joint space `XY`, then marginal visitation are
    obtained from the joint bin visits.
- [`SymbolicPermutation`](@ref). Each timeseries is separately [`encode`](@ref)d according
    to its ordinal pattern.
- [`Dispersion`](@ref). Each timeseries is separately [`encode`](@ref)d according to its
    dispersion pattern.

Many more implementations are possible. Each new implementation gives one new
way of estimating the [`ContingencyMatrix`](@ref)
"""
function marginal_encodings end

marginal_encodings(est::CountOccurrences, x::AbstractVector...) = x

function marginal_encodings(est::SymbolicPermutation{m}, x::AbstractVector...) where {m}
    return [encode.(Ref(est.encoding), embed(xᵢ, m, est.τ).data) for xᵢ in x]
end

function marginal_encodings(est::Dispersion, x::AbstractVector...)
    return [symbolize_for_dispersion(est, xᵢ) for xᵢ in x]
end

function marginal_encodings(
        est::ValueHistogram{<:FixedRectangularBinning{D}},
        xs::AbstractVector...) where D

    ϵmin = est.binning.ϵmin[1]
    ϵmax = est.binning.ϵmax[1]
    N = est.binning.N

    encoders = [RectangularBinEncoding(FixedRectangularBinning(ϵmin, ϵmax, N, 1), x) for x in xs]
    [encode.(Ref(e), x) for (x, e) in zip(xs, encoders)]
end
