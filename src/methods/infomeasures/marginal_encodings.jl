import ComplexityMeasures: symbolize_for_dispersion
export marginal_encodings

"""
    marginal_encodings(est::ProbabilitiesEstimator, x::VectorOrDataset...)

Encode/discretize each input vector `xᵢ ∈ x` according to a procedure determined by `est`.
Any `xᵢ ∈ X` that are multidimensional ([`Dataset`](@ref)s) will be encoded column-wise,
i.e. each column of `xᵢ` is treated as a timeseries and is encoded separately.

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

function marginal_encodings(est, x::VectorOrDataset...)
    return marginally_encode_variable.(Ref(est), x)
end

function marginally_encode_variable(est, x::AbstractDataset)
    return Dataset(marginally_encode_variable.(Ref(est), columns(x))...)
end

function marginally_encode_variable(est::CountOccurrences, x::AbstractVector)
    return x
end

function marginally_encode_variable(est::SymbolicPermutation{m}, x::AbstractVector) where {m}
    emb = embed(x, m, est.τ).data
    return encode.(Ref(est.encoding), emb)
end

function marginally_encode_variable(est::Dispersion, x::AbstractVector)
    return symbolize_for_dispersion(est, x)
end

function marginally_encode_variable(
        est::ValueHistogram{<:FixedRectangularBinning{D}},
        x::AbstractVector) where D
    ϵmin = est.binning.ϵmin[1]
    ϵmax = est.binning.ϵmax[1]
    N = est.binning.N
    encoder = RectangularBinEncoding(FixedRectangularBinning(ϵmin, ϵmax, N, 1))
    return encode.(Ref(encoder), x)
end
