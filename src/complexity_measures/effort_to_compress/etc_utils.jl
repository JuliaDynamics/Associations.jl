using StatsBase
using Entropies
export get_windows

# TODO: fix indexing and re-initialization, so this function can be re-used without 
# providing newly initialized `frequency_of_symbols`, `unique_symbols` and `sorted_sequence`
# It works for now, but initialization of the arrays need to happen *before* calling this 
# function.
"""
estimate_histogram!(frequency_of_symbols, unique_symbols, sorted_sequence, kmax::Int) → Int

Count the occurrences of each of the unique symbols in `sorted_sequence`. The results
are written into the pre-allocated arrays `frequency_of_symbols` and `unique_symbols`.
Assummes `length(sorted_sequence) >= 2`.

This is used for determining which symbols to substitute in a given iteration of the 
[`effort_to_compress`](@ref) algorithm. Because the [`effort_to_compress`](@ref) algorithm 
is intented to be used repeatedly in real applications, pre-allocation is necessary for 
efficiency. Therefore, this function assumes a fixed-length pre-allocated 
`sorted_sequence`-vector, for which we consider the elements at indices `1:kmax` (as 
symbols are substituted in [`effort_to_compress`](@ref), `kmax` decreases).

Returns the number of unique elements found in `sorted_sequence`.
"""
function estimate_histogram!(frequency_of_symbols, unique_symbols, sorted_sequence, 
    kmax::Int)
    frequency_of_symbols .= 0
    n_unique = 1
    
    unique_symbols[n_unique] = sorted_sequence[n_unique]
    frequency_of_symbols[n_unique] = 1

    i = 1
    while i < kmax
        i += 1
        last = sorted_sequence[i - 1]
        next = sorted_sequence[i]
        if last == next
            frequency_of_symbols[n_unique] += 1
        else
            n_unique += 1
            unique_symbols[n_unique] = next
            frequency_of_symbols[n_unique] += 1
        end
    end

    #@show "histogram: ", frequency_of_symbols, unique_symbols, sorted_sequence

end


"""
    ReplacementStrategyInCaseOfEquality

How to decide which symbol to replace in case two or more symbols have the same frequency.
"""
abstract type ReplacementStrategyInCaseOfEquality end

""" 
    PickFirst <: ReplacementStrategyInCaseOfEquality

Pick the first element (according to the sorting algorithm used).
"""
struct PickFirst <: ReplacementStrategyInCaseOfEquality end



"""
Return the sequence where the element with the highest frequency is replaced by 
the new symbol `replacement_symbol`.
"""
function substituted_sequence(x::AbstractVector{T}, 
        symbol_sequence, frequency_of_symbols, unique_symbols, 
        kmax::Int, replacement_strategy::PickFirst) where T <: Integer
    
    replacement_value = maximum(x) + 1

    # The new x will always be smaller than the original x.
    # We don't know how much smaller, so reserve sufficient space.
    new_x = Vector{T}(undef, 0)
    sizehint!(new_x, length(x))
    
    idx_max = argmax(@view frequency_of_symbols[1:kmax])
    symbol_to_replace = (@view unique_symbols[1:kmax])[idx_max]
    indices_to_replace = @views findall(xᵢ -> xᵢ == symbol_to_replace, symbol_sequence[1:kmax])
    j = 1

    while j <= length(x)
        if j in indices_to_replace
            push!(new_x, replacement_value)
            j += 2
        else
            if (j <= length(x))
                push!(new_x, x[j])
                j += 1
            end
        end
    end

    # Shrink new vector to its size
    sizehint!(new_x, length(new_x))

    return new_x
end

"""
Fill the pre-allocated `symbol_sequence` vector with two-letter symbols made from the 
elements `x[1:kmax+1]`. Also, sort the sequence and write in-place to the pre-allocated
`sorted_sequence` vector.
"""
function make_symbols!(symbol_sequence, sorted_sequence, x, kmax::Int) 
    for j = 1:kmax
        symbol_sequence[j] = make_symbol(x, j)
    end


    sorted_sequence[1:kmax] .= symbol_sequence[1:kmax]
    sort!(sorted_sequence, alg = Base.Sort.InsertionSort)
end

make_symbol(x, i::Int) = @views "$(x[i])$(x[i+1])"
make_symbol_svec(x, i::Int) = @views SVector(x[i], x[i+1])

function get_windows(x, window_length, step::Int = 1)
    [i:i+window_length for i in 1:step:length(x)-window_length]
end