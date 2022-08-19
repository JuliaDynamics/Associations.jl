using Entropies, StatsBase
export get_windows

function compress!(x::AbstractVector{T}, 
        symbol_sequence::Vector{String}, 
        sorted_sequence::Vector{String}, 
        unique_symbols::Vector{String}, 
        frequency_of_symbols::Vector{Int}; 
        replacement_strategy::ReplacementStrategyInCaseOfEquality = PickFirst()
        ) where {T <: Integer}

    if (length(x) == 1)
        return x
    end

    # With two-letter symbols, the last entry of `x` enters the `length(x)-1`-th symbol.
    kmax = length(x) - 1
    make_symbols!(symbol_sequence, sorted_sequence, x, kmax)
    estimate_histogram!(frequency_of_symbols, unique_symbols, sorted_sequence, kmax)

    return substituted_sequence(x, 
        symbol_sequence, 
        frequency_of_symbols, 
        unique_symbols, 
        kmax, 
        replacement_strategy)
end

function compress(x::AbstractVector{T}; 
        replacement_strategy::ReplacementStrategyInCaseOfEquality = PickFirst()
        ) where {T <: Integer}

    if (length(x) == 1)
        return x
    end

    # With two-letter symbols, the last entry of `x` enters the `length(x)-1`-th symbol.
    kmax = length(x) - 1

    symbol_sequence = Vector{String}(undef, kmax)
    sorted_sequence = Vector{String}(undef, kmax)
    unique_symbols = Vector{String}(undef, kmax)
    frequency_of_symbols = zeros(Int, kmax)

    return compress!(x, symbol_sequence, sorted_sequence, unique_symbols, frequency_of_symbols; replacement_strategy=replacement_strategy)
end


function compression_complexity(x::AbstractVector, algorithm::EffortToCompress)
    # With two-letter symbols, the last entry of `x` enters the `length(x)-1`-th symbol.
    kmax = length(x) - 1

    # Pre-allocate for repeated use (in the worst case, `kmax - 1` compressions occur)
    # TODO: pre-allocation doesn't really do anything at the moment. Need to fix 
    # indexing and re-initialization in `compress!` and its sub-functions first.
    # Initializing here for every call will do for now.
    symbol_sequence = Vector{String}(undef, kmax)
    sorted_sequence = Vector{String}(undef, kmax)
    unique_symbols = Vector{String}(undef, kmax)
    frequency_of_symbols = zeros(Int, kmax)

    # The ETC complexity measure records how many compression steps it takes to reach
    # a zero-entropy (a single-element or constant) symbol sequence. We use a float here
    # to ensure type-stability if the value is to be normalized.
    N = 0.0

    # Keep track of how many of the sequence symbols have been compressed.
    n_reduced = 0

    new_x = copy(x)
    while genentropy(Dataset(new_x), CountOccurrences()) > 0 && length(new_x) > 1
        symbol_sequence .= ""
        sorted_sequence .= ""
        unique_symbols .= ""
        frequency_of_symbols .= 0
        new_x = compress(new_x)
        N += 1
    end

    return algorithm.normalize ? N/kmax : N
end

function compression_complexity(x::AbstractVector{T}, 
        algorithm::EffortToCompressSlidingWindow) where {T <: Integer}
    windows = get_windows(x, algorithm.window_size, algorithm.step)
    alg = EffortToCompress(normalize = algorithm.normalize)
    return @views [compression_complexity(x[window], alg) for window in windows]
end
