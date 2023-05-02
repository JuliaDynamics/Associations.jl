using StatsBase: countmap
using StateSpaceSets: AbstractStateSpaceSet
using StaticArrays: SVector, SA
export get_windows


allidentical(x) = all(xᵢ -> xᵢ == first(x), x)

haszeroentropy(x) = allidentical(x) || length(x) == 1

histogram(x) = countmap(x)

function vecints_to_int(x::AbstractVector{J}, base) where J <: Integer
    s = zero(J)
    for xᵢ in x
       s = s * base + xᵢ
    end
    return s
end

"""
    sequential_symbol_pairs(x::AbstractVector{J}) where J <: Integer → Vector{SVector{2, J}}

Create a sequence of two-letter symbol pairs from an integer valued time series.
"""
function sequential_symbol_pairs(x::AbstractVector{J}) where {J}
    pairs = zeros(SVector{2, J}, length(x) - 1)

    for i in 1:length(x) - 1
        pairs[i] = SA[x[i], x[i+1]]
    end
    return pairs
end

"""
    sequential_symbol_pairs(x::AbstractVector{J}, y::AbstractVector{J}) where {J <: Integer} → Vector{SVector{2, Tuple{J, J}}}

Generate symbol pairs from the integer-valued, same-length vectors `x` and `y`.
The `i`-th element of the returned vector is the `SVector`
`SA[(x[i], x[i+1]), (y[i], y[i+1])]`.
"""
function sequential_symbol_pairs(x::AbstractVector{J}, y::AbstractVector{J}) where {J <: Integer}
    length(x) == length(y) || throw(ArgumentError("lengths of `x` and `y` must be equal"))
    length(x) >= 2 || throw(ArgumentError("Lengths of x and y must be >= 2"))

    joint_pairs = Vector{SVector{2, Tuple{Int64, Int64}}}(undef, 0)
    for i = 1:length(x) - 1
        push!(joint_pairs, SA[(x[i], x[i+1]), (y[i], y[i+1])])
    end
    return joint_pairs
end

"""
    symbol_sequence(x::AbstractStateSpaceSet{D, T}, alphabet_size) where {D <: Integer, T <: Integer} ⇥ Vector{Int}

Create an encoded one-dimensional integer symbol sequence from a integer-valued , assuming that the alphabet used for the original symbolization
has `alphabet_size` possible symbols.

In the context of the effort-to-compress algorithm, `x` would typically
be a sub-StateSpaceSet as obtained by `view(d::StateSpaceSet, indices)`.
"""
function symbol_sequence(x::AbstractStateSpaceSet{D, T}, alphabet_size::J) where {D, T, J <: Integer}
    return [vecints_to_int(xᵢ, alphabet_size) for xᵢ in x]
end

# function get_windows(x, window_length::Int, step::Int = 1)
#     [i:i+window_length for i in 1:step:(length(x) - window_length)]
# end

# get_windows(x, sw::ConstantWidthSlidingWindow{<:CompressionComplexityEstimator}) = 
#     get_windows(x, sw.width, sw.step)
