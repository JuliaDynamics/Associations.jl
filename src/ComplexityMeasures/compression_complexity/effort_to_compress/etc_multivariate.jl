using DelayEmbeddings
export symbol_sequence


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
