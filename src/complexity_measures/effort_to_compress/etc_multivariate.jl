using DelayEmbeddings
export symbol_sequence

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
function symbol_sequence(x::AbstractDataset{D, T}, alphabet_size) where {D, T}
    encoded_sequence = [vecints_to_int(xᵢ, alphabet_size) for xᵢ in x]
end


function compression_complexity(x::AbstractDataset{D, T}, algorithm::EffortToCompress) where {D, T}
    !isnothing(algorithm.alphabet_size) || throw(ArgumentError("Alphabet size must be specified when input is multivariate"))
    algorithm.alphabet_size >= 2 || throw(ArgumentError("Alphabet size must be at least 2."))
    encoded_sequence = symbol_sequence(x, algorithm.alphabet_size)
    return compression_complexity(encoded_sequence, algorithm)
end

function compression_complexity(x::AbstractDataset{D, T}, algorithm::EffortToCompressSlidingWindow) where {D, T}
    !isnothing(algorithm.alphabet_size) || throw(ArgumentError("Alphabet size must be specified when input is multivariate"))
    algorithm.alphabet_size >= 2 || throw(ArgumentError("Alphabet size must be at least 2."))
    
    encoded_sequence = symbol_sequence(x, algorithm.alphabet_size)
    windows = get_windows(x, algorithm.window_size, algorithm.step)
    alg = EffortToCompress(normalize = algorithm.normalize)
    
    return @views [compression_complexity(encoded_sequence[w], alg) for w in windows]
end
