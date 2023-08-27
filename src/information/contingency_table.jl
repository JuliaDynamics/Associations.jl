import ComplexityMeasures: probabilities
using ComplexityMeasures: RelativeAmount
using NamedArrays
import NamedArrays: Name

export contingency_table
export ContingencyTable
export Name

# So that we can mix discrete-valued state space sets with discrete-valued regular
# vectors.
unique_elements(x) = unique(x)
unique_elements(x::AbstractStateSpaceSet) = unique(x.data)

"""
    ContingencyTable(cts, outcomes)

A contingency table that consists of the counts `cts` (an `Array{N, Int}`, where `N`
is the dimension) and the `outcomes` (an `N`-tuple), where `outcomes[i]` is the
outcomes for the `i`-th dimension.

To actually construct a contingency table, use the [`contingency_table`](@ref) function.

## Indexing

Values of a `ContingencyTable` can be accessed using traditional integer-based indexing.
However, since the counts are internally implemented as a `NamedArray`, you can
also use the *outcomes* to index the contingency table. If so, you must wrap the
value with which you're indexing with `Name` (see example below).


## Example

```julia
using StateSpaceSets
using CausalityTools
using Random; rng = MersenneTwister(1234)

n = 100
x = rand(["a", "b", "c", 2], n)
y = StateSpaceSet(rand(rng, n, 2))
encodings = PerPointEncoding(CategoricalEncoding(x), OrdinalPatternEncoding(3))
c = contingency_table(encodings, x, y)

z = StateSpaceSet(randn(n, 2))
encodings = PerPointEncoding(
    CategoricalEncoding(x),
    OrdinalPatternEncoding(2),
    GaussianCDFEncoding{2}(μ = 0.0, σ = 0.3; c = 3)
)
c = contingency_table(encodings, x, y, z)

# Subset a two-dimensional marginal using traditional indexing
c[:, :, 3]

# The same, but using the outcome as the index
c[:, :, Name(c.outcomes[3][3])]

# Subsetting a one-dimensional marginal.
c[:, Name(SVector(2, 1)), Name(c.outcomes[3][3])]
```
"""
struct ContingencyTable{T, N, V} <: AbstractArray{T, N}
    cts::NamedArray{T, N}
    outcomes::V
end

Base.show(io::IO, c::ContingencyTable) = show(io, summary_string(c))
function Base.show(io::IO, m::MIME"text/plain", c::ContingencyTable)
    arr = repr("text/plain", c.cts)
    # remove first two lines (these just show `NamedArray` info)
    arr_final = join(split(arr, '\n')[3:end], "\n  ")
    total_str = join([summary_string(c), arr_final], '\n')
    print(io,  total_str)
end

function summary_string(c::ContingencyTable)
    ds = size(c.cts)
    n = sum(c.cts)
    summary_string = join(["$d" for d in ds], 'x') * " ContingencyTable with n=$(n) counts"
    return summary_string
end

Base.getindex(c::ContingencyTable, i...) = getindex(c.cts, i...)
Base.setindex(c::ContingencyTable, i...) = setindex(c.cts, i...)
Base.iterate(c::ContingencyTable, state=1) = iterate(c.cts, state)
Base.eltype(c::ContingencyTable) = eltype(c.cts)
Base.length(c::ContingencyTable) = length(c.cts)
Base.size(c::ContingencyTable) = size(c.cts)
Base.firstindex(c::ContingencyTable) = firstindex(c.cts)
Base.lastindex(c::ContingencyTable) = lastindex(c.cts)

"""
    contingency_table(x₁, x₂, ..., xₙ) → ContingencyTable{N}

Construct an `N`-dimensional contingency table of counts from the input vectors
``x_1, x_2, \\ldots, x_N``, where each ``xₖ` must be an iterable containing
discrete values.

These discrete iterables are typically `Vector{Int}` constructed from input data using
[`encode`](@ref) in combination with some [`Discretization`](@ref).
"""
function contingency_table(x...)

    # Get marginal probabilities and outcomes
    L = length(x)
    cts, lmaps, encoded_outcomes = counts_table(x...)
    # lmaps[i]: a `Dict{outcome_type, Int}` containing the conversion between the
    #   internally encoded outcomes for the `i`-th input, and the actual outcomes
    #   for the `i`-th input.
    actual_outcomes = map(i -> to_outcomes(lmaps[i], encoded_outcomes[i]), tuple(1:L...))
    cts_named = NamedArray(cts, actual_outcomes)
    return ContingencyTable(cts_named, actual_outcomes)
end

function to_outcomes(lmap::Dict, encoded_outcomes::Vector{<:Integer})
    # We want the encoded integers as keys and the actual outcomes as values.
    lmap_swapped = Dict(values(lmap) .=> keys(lmap))
    return [lmap_swapped[ωᵢ] for ωᵢ in encoded_outcomes]
end

# For this to work generically, we must map unique elements to integers.
function contingency_matrix(x...)
    # TODO: Inverse map from integer-encoded outcomes to the original outcomes.
    # marginal_outcomes = [map(k -> lmap[k], last(pΩ)) for (pΩ, lmap) in zip(pΩs, lmaps)]
    cts, lmaps, outcomes
    probs = cts ./ L # relative frequency estimation.
    return ContingencyMatrix(
        probs,
        cts,
    )
end

function counts_table(x...)
    Ls = length.(x);
    if !allequal(Ls)
        throw(ArgumentError("Input data must have equal lengths. Got lengths $Ls."))
    end
    L = first(Ls)

    # Map the input data to integers. This ensures compatibility with *any* input type.
    # Then, we can simply create a joint `StateSpaceSet{length(x), Int}` and use its elements
    # as `CartesianIndex`es to update counts.
    lvl = tolevels.(x)
    levels = (first(l) for l in lvl)
    lmaps = [last(l) for l in lvl]

    # Create the table with correct dimensions, assumming the outcome space is
    # fully determined by the elements that are present in `x`.
    table_dims = length.(unique_elements.(x));
    cts = zeros(Int, table_dims)
    # Each element in `X` isa `SVector{m, Int}`, so can be treated as a cartesian index.
    X = StateSpaceSet(levels...)

    # We sort, so that the positions in `cts` will correspond to the indices on
    # each of the axes of `cts`. Note: these are not the *actual* outcomes, but the
    # internal integer representation of each outcome. We need to use `lmaps` to convert
    # back.
    for ix in X
        cts[ix...] += 1
    end

    # One set of outcomes per input
    outcomes = sort!.(unique!.(columns(X)))
    return cts, lmaps, outcomes
end

function to_cartesian(x)
    (CartesianIndex.(xᵢ...) for xᵢ in x)
end

"""
    tolevels!(levels, x) → levels, dict
    tolevels(x) → levels, dict

Apply the bijective map ``f : \\mathcal{Q} \\to \\mathbb{N}^+`` to each `x[i]` and store
the result in `levels[i]`, where `levels` is a pre-allocated integer vector such that
`length(x) == length(levels)`.

``\\mathcal{Q}`` can be any space, and each ``q \\in \\mathcal{Q}`` is mapped to a unique
integer  in the range `1, 2, …, length(unique(x))`. This is useful for integer-encoding
categorical data such as strings, or other complex data structures.

The single-argument method allocated a `levels` vector internally.

`dict` gives the inverse mapping.
"""
function tolevels!(levels, x)
    @assert length(levels) == length(x)
    lmap = _levelsmap(x)
    for i in eachindex(x)
        levels[i] = lmap[x[i]]
    end
    return levels, lmap
end

function tolevels(x)
    lmap = _levelsmap(x)
    levels = zeros(Int, length(x))
    for i in eachindex(x)
        levels[i] = lmap[x[i]]
    end
    return levels, lmap
end

# Ugly hack, because levelsmap doesn't work out-of-the-box for statespacesets.
_levelsmap(x) = levelsmap(x)
_levelsmap(x::AbstractStateSpaceSet) = levelsmap(x.data)
