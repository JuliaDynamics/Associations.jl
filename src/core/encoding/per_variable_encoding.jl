import ComplexityMeasures: codify, encode
import ComplexityMeasures: OutcomeSpace

using DelayEmbeddings: embed
export VariableEncoding
export codify

# TODO: implement this generically for `Encodings` too (will require type-parameterized
# number of elements for the Encodings).

"""
    VariableEncoding <: Discretization
    VariableEncoding(encoding::Encoding)

Given multiple dataset `xs::AbstractStateSpaceSet`, [`encode`](@ref) the `i`-th
variable/column using the given `encoding`. The same encoding is applied to each column to
ensure that each column has the same length after encoding.

Encoding is done internally using [`codify`](@ref), which typically runs a sliding
window (of width dictated by the given `encoding`) across each variable, encoding
each window to a unique integer.
"""
struct VariableEncoding{N} <: Discretization{N}
    outcome_spaces::NTuple{N, OutcomeSpace}
    function VariableEncoding(outcome_spaces::NTuple{N, OutcomeSpace}) where N
        if N > 1
            s = "It is currently only possible to use the same `OutcomeSpace` for all " *
                "variables. Got $N different encodings"
            throw(ArgumentError(s))
        end
        new{N}(outcome_spaces)
    end
end
# TODO: we can use encoded_space_cardinality if multiple outcome spaces are used
# to test if the outcome spaces can be matched.

function VariableEncoding(o::OutcomeSpace)
    return VariableEncoding((o,))
end

function encode(encoding::VariableEncoding{1}, x::Vararg{Any, 1})
    e = first(encoding.outcome_spaces)
    x̂ = codify(e, first(x))
    return x̂::Vector{<:Integer}
end

function encode(encoding::VariableEncoding{1}, x::Tuple)
    e = first(encoding.outcome_spaces)
    x̂ = map(xᵢ -> codify(e, xᵢ), x)
    N = length(x)
    return x̂::NTuple{N, Vector{<:Integer}}
end

function encode(encoding::VariableEncoding{1}, x::AbstractStateSpaceSet)
    return encode(encoding, columns(x))
end
