import ComplexityMeasures: symbolize_for_dispersion
export marginal_encodings

"""
    marginal_encodings(est::ProbabilitiesEstimator, x::VectorOrStateSpaceSet...)

Encode/discretize each input vector `xᵢ ∈ x` according to a procedure determined by `est`.
Any `xᵢ ∈ X` that are multidimensional ([`StateSpaceSet`](@ref)s) will be encoded column-wise,
i.e. each column of `xᵢ` is treated as a timeseries and is encoded separately.

This is useful for computing any discrete information theoretic quantity, and is
used internally by [`contingency_matrix`](@ref).

## Supported estimators

- [`ValueHistogram`](@ref). Bin visitation frequencies are counted in the joint space `XY`,
    then marginal visitations are obtained from the joint bin visits.
    This behaviour is the same for both [`FixedRectangularBinning`](@ref) and
    [`RectangularBinning`](@ref) (which adapts the grid to the data).
    When using [`FixedRectangularBinning`](@ref), the range along the first dimension
    is used as a template for all other dimensions.
- [`OrdinalPatterns`](@ref). Each timeseries is separately [`encode`](@ref)d according
    to its ordinal pattern.
- [`Dispersion`](@ref). Each timeseries is separately [`encode`](@ref)d according to its
    dispersion pattern.

Many more implementations are possible. Each new implementation gives one new
way of estimating the [`ContingencyMatrix`](@ref)
"""
function marginal_encodings end

function marginal_encodings(est, x::VectorOrStateSpaceSet...)
    return marginally_encode_variable.(Ref(est), x)
end

function marginally_encode_variable(est, x::AbstractStateSpaceSet)
    return StateSpaceSet(marginally_encode_variable.(Ref(est), columns(x))...)
end

function marginally_encode_variable(est::CountOccurrences, x::AbstractVector)
    return x
end

function marginally_encode_variable(est::OrdinalPatterns{m}, x::AbstractVector) where {m}
    emb = embed(x, m, est.τ).data
    return encode.(Ref(est.encoding), emb)
end

function marginally_encode_variable(est::Dispersion, x::AbstractVector)
    return symbolize_for_dispersion(est, x)
end

function marginally_encode_variable(
        est::ValueHistogram{<:FixedRectangularBinning{D}},
        x::AbstractVector) where D
    range = first(est.binning.ranges)
    ϵmin = minimum(range)
    ϵmax = maximum(range)
    N = length(range)
    encoder = RectangularBinEncoding(FixedRectangularBinning(ϵmin, ϵmax, N, 1))
    return encode.(Ref(encoder), x)
end

# Special treatment for RectangularBinning. We create the joint embedding, then
# extract marginals from that. This could probably be faster,
# but it *works*. I'd rather things be a bit slower than having marginals
# that are not derived from the same joint distribution, which would hugely increase
# bias, because we're not guaranteed cancellation between entropy terms
# in higher-level methods.
function marginal_encodings(est::ValueHistogram{<:RectangularBinning}, x::VectorOrStateSpaceSet...)
    X = StateSpaceSet(StateSpaceSet.(x)...)
    encoder = RectangularBinEncoding(est.binning, X)

    bins = [vec(encode_as_tuple(encoder, pt))' for pt in X]
    joint_bins = reduce(vcat, bins)
    idxs = size.(x, 2) #each input can have different dimensions
    s = 1
    encodings = Vector{Vector}(undef, 0)
    for (i, cidx) in enumerate(idxs)
        variable_subset = s:(s + cidx - 1)
        s += cidx
        y = @views joint_bins[:, variable_subset]
        for j in size(y, 2)
            push!(encodings, y[:, j])
        end
    end

    return encodings
end

# A version of `cartesian_bin_index` that directly returns the joint bin encoding
# instead of converting it to a cartesian index.
function encode_as_tuple(e::RectangularBinEncoding, point::SVector{D, T}) where {D, T}
    ranges = e.ranges
    if e.precise
        # Don't know how to make this faster unfurtunately...
        bin = map(searchsortedlast, ranges, point)
    else
        bin = floor.(Int, (point .- e.mini) ./ e.widths) .+ 1
    end
    return bin
end
