"""
    lexicographical_subsets(set::Vector{Tuple{Int, Int}}, cardinality::Int; 
        exclude::Tuple{Int, Int} = nothing)

Generating lexicographical subsets with `cardinality` from the given `set` of 
integer tuples of the form `(i, τ)`, where `i` is the index of the variable
(in the original `StateSpaceSet` given to `infer_graph`) , and `τ` is the 
lag. If a tuple is given as `exclude`, then this `(index, lag)` is excluded
from the vector of returned subsets.
"""
function lexicographical_subsets(set::Vector{Tuple{Int,Int}}, cardinality::Int; exclude::Tuple{Int,Int}=nothing)
    filtered_set = isnothing(exclude) ? set : filter(x -> x != exclude, set)
    return [subset for subset in combinations(filtered_set, cardinality)]
end


"""
    fdr_adjust(alg::PCMCI, p_values::SortedDict{Tuple{Int, Int, Int}, Float64})

Adjust `p_values` to control false discovery rate using the Benjamini-Hochberg (1995)
procedure. 

## Arguments
- **`p_values`*:  Dictionary of p-values with keys representing (i, j, τ)
- **`α`**: Significance level.
"""
function fdr_adjust(alg::PCMCI, p_values::SortedDict{Tuple{Int,Int,Int},Float64})
    # Extract p-values and sort them in ascending order
    sorted_keys = collect(keys(p_values))
    sorted_p_values = collect(values(p_values))
    n_tests = length(sorted_p_values)

    # Sort the p-values and keep track of the original indices
    sorted_indices = sortperm(sorted_p_values)
    sorted_p_values = sorted_p_values[sorted_indices]
    sorted_keys = sorted_keys[sorted_indices]

    # Apply the BH procedure to compute adjusted p-values
    adjusted_p_values = Vector{Float64}(undef, n_tests)
    for i in 1:n_tests
        rank = i
        # Compute the BH critical value
        adjusted_p_values[i] = min(1.0, sorted_p_values[i] * n_tests / rank)
    end

    # Ensure adjusted p-values are non-decreasing
    for i in (n_tests-1):-1:1
        adjusted_p_values[i] = min(adjusted_p_values[i], adjusted_p_values[i+1])
    end

    # Return a dictionary with similar structure, but with the adjusted p-values
    adjusted_p_values_dict = SortedDict{Tuple{Int,Int,Int},Float64}()
    for i in 1:n_tests
        adjusted_p_values_dict[sorted_keys[i]] = adjusted_p_values[i]
    end

    return adjusted_p_values_dict
end