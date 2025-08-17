using ComplexityMeasures
export codified_marginals

"""
    codified_marginals(o::OutcomeSpace, x::VectorOrStateSpaceSet...)

Encode/discretize each input vector (e.g. timeseries) `xᵢ ∈ x` according to a procedure
determined by `o`. 

For some outcome spaces, the encoding is sequential (i.e. time ordering matters). 
Any `xᵢ ∈ X` that are multidimensional ([`StateSpaceSet`](@extref StateSpaceSets.StateSpaceSet)s) will be encoded
column-wise, i.e. each column of `xᵢ` is treated as a timeseries and is encoded separately.

This is useful for discretizing input data when computing some 
[`MultivariateInformationMeasure`](@ref). This method is used internally by
both the [`JointProbabilities`](@ref) and [`EntropyDecomposition`](@ref) estimators
to handle discretization.

## Supported estimators

- [`ValueBinning`](@extref ComplexityMeasures.ValueBinning). Bin visitation frequencies are counted in the joint space `XY`,
    then marginal visitations are obtained from the joint bin visits.
    This behaviour is the same for both [`FixedRectangularBinning`](@extref ComplexityMeasures.FixedRectangularBinning) and
    [`RectangularBinning`](@extref ComplexityMeasures.RectangularBinning) (which adapts the grid to the data).
    When using [`FixedRectangularBinning`](@extref ComplexityMeasures.FixedRectangularBinning), the range along the first dimension
    is used as a template for all other dimensions.
- [`OrdinalPatterns`](@extref ComplexityMeasures.OrdinalPatterns). Each timeseries is separately [`codify`](@ref)-ed by 
    embedding the timeseries, then sequentially encoding the ordinal patterns of 
    the embedding vectors.
- [`Dispersion`](@extref ComplexityMeasures.Dispersion). Each timeseries is separately [`codify`](@ref)-ed by 
    embedding the timeseries, then sequentially encoding the embedding vectors
    according to their dispersion pattern (which for each embedding vector is computed
    relative to all other embedding vectors).
- [`CosineSimilarityBinning`](@extref ComplexityMeasures.CosineSimilarityBinning). Each timeseries is separately [`codify`](@ref)-ed
    by embedding the timeseries, the encoding the embedding points in a 
    in a sequential manner according to the cosine similarity of the embedding vectors.
- [`UniqueElements`](@extref ComplexityMeasures.UniqueElements). Each timeseries is [`codify`](@ref)-ed according to 
    its unique values (i.e. each unique element gets assigned a specific integer).

More implementations are possible.
"""
function codified_marginals end

function codified_marginals(d::CodifyVariables, x::VectorOrStateSpaceSet...)
    T = eltype(d.outcome_spaces) # assume identical outcome spaces.
    if !allequal(typeof.(d.outcome_spaces))
        throw(ArgumentError("Outcome space for each marginal must be identical. Got outcome spaces of type $T"))
    end
    o = first(d.outcome_spaces) # we can do this because we assume all out come spaces are the same
    return codified_marginals(o, x...)
end

function codified_marginals(o::OutcomeSpace, x::VectorOrStateSpaceSet...)
    return codify_marginal.(Ref(o), x)
end

# Generic dispatch to ComplexityMeasures.jl. We override if something special
# needs to happen. For example, for ValueBinning we override such that 
# we bin in the joint space to reduce bias.
function codify_marginal(o::OutcomeSpace, x::VectorOrStateSpaceSet)
    return codify(o, x)
end
# Apply per column.
function codify_marginal(o::OutcomeSpace, x::AbstractStateSpaceSet)
    return StateSpaceSet(codify_marginal.(Ref(o), columns(x))...)
end

# ------------------------------------------------------------------------
# Outcome space specific implementations
# ------------------------------------------------------------------------

# TODO: maybe construct a convenience wrapper where the user can avoid constructing the
# joint space, for performance benefits (but increased bias).
function codify_marginal(
    o::ValueBinning{<:FixedRectangularBinning{D}},
    x::AbstractVector) where D
    range = first(o.binning.ranges)
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
function codified_marginals(o::ValueBinning{<:RectangularBinning}, x::VectorOrStateSpaceSet...)
    # TODO: The following line can be faster by explicitly writing out loops that create the 
    # joint embedding vectors.
    X = StateSpaceSet(StateSpaceSet.(x)...)
    encoder = RectangularBinEncoding(o.binning, X)

    bins = [vec(encode_as_tuple(encoder, pt))' for pt in X]
    joint_bins = reduce(vcat, bins)
    idxs = size.(x, 2) #each input can have different dimensions
    s = 1
    encodings = Vector{Vector}(undef, 0)
    for (i, cidx) in enumerate(idxs)
        variable_subset = s:(s+cidx-1)
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
function encode_as_tuple(e::RectangularBinEncoding, point::SVector{D,T}) where {D,T}
    ranges = e.ranges
    if e.precise
        # Don't know how to make this faster unfurtunately...
        bin = map(searchsortedlast, ranges, point)
    else
        bin = floor.(Int, (point .- e.mini) ./ e.widths) .+ 1
    end
    return bin
end