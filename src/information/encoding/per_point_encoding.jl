using NamedArrays
using ComplexityMeasures: Encoding
import ComplexityMeasures: encode

export PerPointEncoding
export encode

""" A discretization of `N` input datasets. """
abstract type Discretization{N} end

"""
    PerPointEncoding{N}
    PerPointEncoding(encodings::NTuple{N, Encoding})

Given multiple dataset `xs::StateSpaceSet`, [`encode`](@ref) every `mₖ`-dimensional point
`x[k][i]` as an integer using the given
[`Encoding`](@ref)s `encodings[k]`.
"""
struct PerPointEncoding{N} <: Discretization{N}
    encodings::NTuple{N, Encoding}
    function PerPointEncoding(encodings::NTuple{N, Encoding}) where N
        if !(N ≥ 1)
            throw(ArgumentError("PerPointEncoding requires at least 1 dimensions"))
        end
        new{N}(encodings)
    end
end
Base.getindex(e::PerPointEncoding, i) = getindex(e.encodings, i)

function PerPointEncoding(encodings::Vararg{Encoding, N}) where N
    return PerPointEncoding(tuple(encodings...))
end

"""
    encode(encoding::PerPointEncoding{N}, x::Vararg{<:AbstractStateSpaceSet, N})

Encode

## Examples

```julia
x = StateSpaceSet(rand(100, 2))
y = StateSpaceSet(rand(100, 3))
z = StateSpaceSet(rand(100, 4))

# For `x`, we use a relative mean encoding.
ex = RelativeMeanEncoding(0.0, 1.0, n = 3)
# For `y`, we use a combination encoding.
ey = CombinationEncoding(RelativeMeanEncoding(0.0, 1.0, n = 3), OrdinalPatternEncoding(3))
# For `z`, we use ordinal patterns to encode.
ez = OrdinalPatternEncoding(4)

# Encoding two input datasets gives a 2-tuple of Vector{Int}
encode(PerPointEncoding(ex, ey), x, y)

# Encoding three input datasets gives a 3-tuple of Vector{Int}
encode(PerPointEncoding(ex, ey, ez), x, y, z)
"""
function encode(encoding::PerPointEncoding{N}, x...) where {N}
    M = length(x)

    if N != M
        s = "The given `encoding` is for $N input datasets. $M input datasets were given."
        throw(ArgumentError(s))
    end
    Ls = length.(x)
    if !allequal(Ls)
        throw(ArgumentError("All input datasets must have the same length."))
    end
    L = maximum(Ls)

    # x̂ := Encoded `x`s, where `x̂[k]` is the encoded version of `x[k]`
    x̂ = [encode_individual_dataset(encoding[k], x[k]) for k = 1:M]
    return x̂
end

function encode_individual_dataset(encoding::Encoding, x)
    if !(typeof(x) <: AbstractStateSpaceSet)
        encoding = CategoricalEncoding(x)
        x̂ = encode.(Ref(encoding), x)
        return x̂
    end

    # x̂[i] := the integer code for the state vector `x[i]`.
    x̂ = zeros(Int, length(x))
    @inbounds for i in eachindex(x)
        x̂[i] = encode(encoding, x[i])
    end
    return x̂
end

 # The decoding step on the second-to-last line is not possible without actually providing
 # the encodings. Therefore, we need to override the generic implementation of
 # `contingency_table`.
function contingency_table(encoding::PerPointEncoding, x...)
    # This converts each dataset `x[i]::StateSpaceSet` into `x̂[i]::Vector{Int}`,
    # where `length(x[i]) == length(x̂[i])`.
    x̂ = encode(encoding, x...)
    # lmaps[i]: a `Dict{outcome_type, Int}` containing the conversion between the
    #   internally encoded outcomes for the `i`-th input, and the actual outcomes
    #   for the `i`-th input.
    cts, lmaps, encoded_outcomes = counts_table(x̂...)

    # Actual outcomes (these outcomes map correspond to those in `x̂`).
    # We can't actually decode any further than this.
    L = length(x)
    outcomes = map(i -> to_outcomes(lmaps[i], encoded_outcomes[i]), tuple(1:L...))

    decoded_outcomes = map(i -> decode_outcomes(encoding[i], outcomes[i]), tuple(1:L...))
    cts_named = NamedArray(cts, decoded_outcomes)
    return ContingencyTable(cts_named, decoded_outcomes)
end

function decode_outcomes(encoding::Encoding, outcomes::Vector{<:Integer})
    return decode.(Ref(encoding), outcomes)
end
