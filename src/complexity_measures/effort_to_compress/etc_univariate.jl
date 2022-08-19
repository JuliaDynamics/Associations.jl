using DelayEmbeddings, StaticArrays

function replacement_value(x::Vector{SVector{2, J}}) where {J <: Int}
    max_x = 0
    for xᵢ in x
        if maximum(xᵢ[1]) > max_x
            max_x = maximum(xᵢ[1])
        end
        if maximum(xᵢ[2]) > max_x
            max_x = maximum(xᵢ[2])
        end
    end

    # Which value should it be replaced with?
    replacement_value = max_x + 1
end

# TODO: this can be made much faster by not using dictionaries and using pre-allocation
# and clever indexing. 
"""
    pair_to_be_replaced(x::Vector{SVector{2, J}}) where {J <: Int} → SVector{2, J}

Compute the most frequently occurring symbol pair `xᵢ`, where each `xᵢ ∈ x` is a symbol 
pair. In the case of ties, we return the first pair, according to the underlying (lack of) 
sorting.
"""
function pair_to_be_replaced(x::Vector{SVector{2, J}}) where {J <: Int}
    # The pair to be replaced is the one that occurs most often, so compute occurrences
    # of each unique pair. 
    hist = histogram(x)
    histkeys = keys(hist) |> collect
    frequencies = values(hist) |> collect

    # In the case of frequency ties, we pick the first (the algorithm doesn't care).
    maxfreq_idx = findfirst(x -> x == maximum(frequencies), frequencies)
    pair_to_be_replaced = histkeys[maxfreq_idx]

    return pair_to_be_replaced
end

function non_sequential_recursive_pair_substitution(pairs, pair_to_replace, 
        symbol_to_replace_with)
    
    compressed_x = Vector{Int}(undef, 0)
    sizehint!(compressed_x, length(pairs))

    j = 1
    k = 0

    #println("Total pairs: $(length(pairs))\n$pairs\n---")
    n_replaced = 0
    while j <= length(pairs) && k < 20
        if pairs[j] == pair_to_replace
            # If replacing the j-th symbol pair, then the j-th and (j+1)-th original symbols
            # are merged. Thus, the index must be incremented by 2 (hence "non-sequential").
            push!(compressed_x, symbol_to_replace_with)
            n_replaced += 1

            #print("n_replaced=$n_replaced | Replaced at j=$j ")
            if j == length(pairs) - 1
                j += 1
                @views push!(compressed_x, pairs[j][2])
                j = length(pairs) + 1
            else
                j += 2
            end
            #println("and jumped to j=$j")
        else 
            if j == length(pairs)
                # If at the last symbol pair, keep both symbols.
                @views push!(compressed_x, pairs[j][1])
                @views push!(compressed_x, pairs[j][2])
                #print("n_replaced=$n_replaced | Kept BOTH VALUES OF ORIGINAL PAIR at j=$j ")
                j += 1
                #println(" and jumped to j=$j")
            elseif j < length(pairs)
                    # If at any symbol pair but the last, keep the first symbol.
                    @views push!(compressed_x, pairs[j][1])
                    #print("n_replaced=$n_replaced | Kept original at j=$j ")
                    j += 1
                    #println("and jumped to j=$j")
            end
        end
        
        k += 1
    end
    sizehint!(compressed_x, length(compressed_x))
    return compressed_x
end

function compress(x::Vector{SVector{2, J}}) where {J <: Int}
    # This is a two-letter symbol pair (two integers)
    pair_to_replace = pair_to_be_replaced(x)

    # This is a single integer
    symbol_to_replace_with = replacement_value(x)

    # Compress by using non-sequential recursive pair substitution 
    non_sequential_recursive_pair_substitution(x, pair_to_replace, symbol_to_replace_with)
end



function compression_complexity(x::AbstractVector{J}, algorithm::EffortToCompress) where {J <: Integer}
    seq = copy(x)
    N = 0.0
    while !haszeroentropy(seq)
        seq = compress(symbol_sequence(seq))
        N += 1
    end
    return algorithm.normalize ? N / (length(x) - 1) : N
end


function compression_complexity(x::AbstractVector{T}, 
    algorithm::EffortToCompressSlidingWindow) where {T <: Integer}
    windows = get_windows(x, algorithm.window_size, algorithm.step)
    alg = EffortToCompress(normalize = algorithm.normalize)
    return @views [compression_complexity(x[window], alg) for window in windows]
end

