# This file implements 2-variable joint effort-to-compress (ETC).
# TODO: it is also possible to compute joint ETC for more than two sequences, but it 
# is probably best to use generated functions for this. I'm leaving it for future work.

using StatsBase, Entropies, StaticArrays
function replacement_value(x::Vector{SVector{2, Tuple{J, J}}}, pair_to_replace) where {J <: Integer}
    max_x₁ = zero(J)
    max_x₂ = zero(J)

    for xᵢ in x
        if maximum(xᵢ[1]) > max_x₁
            max_x₁ = maximum(xᵢ[1])
        end
        if maximum(xᵢ[2]) > max_x₂
            max_x₂ = maximum(xᵢ[2])
        end
    end

    return (max_x₁, max_x₂)
end

function pair_to_be_replaced(x::Vector{SVector{2, Tuple{J, J}}}) where {J <: Integer}
    # The pair to be replaced is the one that occurs most often, so compute occurrences
    # of each unique pair. 
    hist = histogram(x)
    histkeys = keys(hist) |> collect
    frequencies = values(hist) |> collect

    # In the case of frequency ties, we pick the first (the algorithm doesn't care).
    maxfreq_idx = findfirst(x -> x == maximum(frequencies), frequencies)

    return histkeys[maxfreq_idx]
end

function non_sequential_recursive_pair_substitution(
        pairs::Vector{SVector{2, Tuple{J, J}}}, 
        pair_to_replace, 
        symbol_to_replace_with) where {J <: Integer}
    #@show pair_to_replace
    #@show symbol_to_replace_with
    T = eltype(first(first(pairs)))
    
    compressed_x = Vector{T}(undef, 0)
    compressed_y = Vector{T}(undef, 0)
    sizehint!(compressed_x, length(pairs))
    sizehint!(compressed_y, length(pairs))

    j = 1
    #println("Total pairs: $(length(pairs))\n$pairs\n---")
    n_replaced = 0
    while j <= length(pairs)
        if pairs[j] == pair_to_replace
            # If replacing the j-th symbol pair, then the j-th and (j+1)-th original symbols
            # are merged. Thus, the index must be incremented by 2 (hence "non-sequential").
            push!(compressed_x, symbol_to_replace_with[1])
            push!(compressed_y, symbol_to_replace_with[2])

            n_replaced += 1

            #print("n_replaced=$n_replaced | Replaced at j=$j ")
            if j == length(pairs) - 1
                j += 1
                @views push!(compressed_x, pairs[j][2][1])
                @views push!(compressed_y, pairs[j][2][2])

                j = length(pairs) + 1
            else
                j += 2
            end
            #println("and jumped to j=$j")
        else 
            if j == length(pairs)
                # If at the last symbol pair, keep both symbols.
                @views push!(compressed_x, pairs[j][1][1])
                @views push!(compressed_x, pairs[j][2][1])
                @views push!(compressed_y, pairs[j][1][2])
                @views push!(compressed_y, pairs[j][2][2])
                #print("n_replaced=$n_replaced | Kept BOTH VALUES OF ORIGINAL PAIR at j=$j ")
                j += 1
                #println(" and jumped to j=$j")
            elseif j < length(pairs)
                # If at any symbol pair but the last, keep the first symbol.
                @views push!(compressed_x, pairs[j][1][1])
                @views push!(compressed_y, pairs[j][1][2])
                #print("n_replaced=$n_replaced | Kept original at j=$j ")
                j += 1
                #println("and jumped to j=$j")
            end
        end
        
    end
    sizehint!(compressed_x, length(compressed_x))
    sizehint!(compressed_y, length(compressed_y))

    return compressed_x, compressed_y
end

function compress(x::AbstractVector{J}, y::AbstractVector{J}) where {J <: Int}
    # This is a sequence of SVector{2, Tuple{Int, Int}}
    seq = symbol_pairs(x, y)

    # This is a two-letter symbol pair (two integers)
    pair_to_replace = pair_to_be_replaced(seq)

    # This is a single integer
    symbol_to_replace_with = replacement_value(seq, pair_to_replace)

    # Compress by using non-sequential recursive pair substitution 
    non_sequential_recursive_pair_substitution(seq, pair_to_replace, symbol_to_replace_with)
end

function compression_complexity(
        x::AbstractVector{J}, 
        y::AbstractVector{J}, 
        algorithm::EffortToCompress
    ) where J <: Integer
    length(x) == length(y) || throw(ArgumentError("lengths of `x` and `y` must be equal"))
    
    # Edge case: one-element vectors return zero regardless of normalization (avoids 
    # division by zero).
    if length(x) == length(y) == 1
        return 0.0
    end

    N = 0.0
    while !haszeroentropy(Dataset(x, y))
        x, y = compress(x, y)
        N += 1
    end

    return algorithm.normalize ? N / (length(x) - 1) : N
end

function compression_complexity(x::AbstractVector{J}, y::AbstractVector{J},
    algorithm::EffortToCompressSlidingWindow) where {J <: Integer}
    windows = get_windows(x, algorithm.window_size, algorithm.step)
    alg = EffortToCompress(normalize = algorithm.normalize)
    return @views [compression_complexity(x[window], y[window], alg) for window in windows]
end
