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

- [`ValueHistogram`](@ref). Bin visitation frequencies are counted in the joint space `XY`,
    then marginal visitations are obtained from the joint bin visits.
    This behaviour is the same for both [`FixedRectangularBinning`](@ref) and
    [`RectangularBinning`](@ref) (which adapts the grid to the data).
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

# Special treatment for RectangularBinning. We create the joint embedding, then
# extract marginals from that. This could probably be faster,
# but it *works*. I'd rather things be a bit slower than having marginals
# that are not derived from the same joint distribution, which would hugely increase
# bias, because we're not guaranteed cancellation between entropy terms
# in higher-level methods.
function marginal_encodings(est::ValueHistogram{<:RectangularBinning}, x::VectorOrDataset...)
    X = Dataset(Dataset.(x)...)
    encoder = RectangularBinEncoding(est.binning, X)
    bins = [vec(encode_as_tuple(encoder, pt))' for pt in X]
    joint_bins = reduce(vcat, bins)

    idxs = size.(x, 2) #each input can have different dimensions
    s = 1
    encodings = Vector{Dataset}(undef, length(idxs))
    for (i, cidx) in enumerate(idxs)
        variable_subset = s:(s + cidx - 1)
        s += cidx
        y = @views joint_bins[:, variable_subset]
        encodings[i] = Dataset(y)
    end

    return encodings

end

# A version of `ComplexityMeasure.encode` that directly returns the joint bin encoding.
function encode_as_tuple(e::RectangularBinEncoding, point)
    (; mini, edgelengths) = e
    # Map a data point to its bin edge (plus one because indexing starts from 1)
    bin = floor.(Int, (point .- mini) ./ edgelengths) .+ 1
    return bin # returns
end
