using StatsBase
using Entropies
export get_windows


allidentical(x) = all(xᵢ -> xᵢ == first(x), x)

haszeroentropy(x) = allidentical(x) || length(x) == 1

histogram(x) = StatsBase.countmap(x)

function vecints_to_int(x::AbstractVector{J}, base) where J <: Integer
    s = zero(J)
    for xᵢ in x
       s = s * base + xᵢ
    end
    return s
end


"""
    symbol_sequence(x::Dataset{D, T}, alphabet_size) where {D <: Integer, T <: Integer}

Create an encoded one-dimensional integer symbol sequence from a `D`-dimensional 
symbol sequence `x`, assuming that the alphabet used for the original symbolization 
has `alphabet_size` possible symbols.

In the context of the effort-to-compress algorithm, `x` would typically 
be a sub-dataset as obtained by `view(d::Dataset, indices)`.
"""
function symbol_sequence(x::AbstractVector{J}) where {J}
    encoded_sequence = zeros(SVector{2, J}, length(x) - 1)

    for i in 1:length(x) - 1
        encoded_sequence[i] = SA[x[i], x[i+1]]
    end
    return encoded_sequence

end


function get_windows(x, window_length, step::Int = 1)
    [i:i+window_length for i in 1:step:length(x)-window_length]
end
