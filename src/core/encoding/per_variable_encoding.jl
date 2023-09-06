using ComplexityMeasures
import ComplexityMeasures: symbolize

using DelayEmbeddings: embed
export PerVariableEncoding
export symbolize

# TODO: implement this generically for `Encodings` too (will require type-parameterized
# number of elements for the Encodings).

"""
    PerVariableEncoding <: Discretization
    PerVariableEncoding(encoding::Encoding)

Given multiple dataset `xs::StateSpaceSet`, [`encode`](@ref) the `i`-th variable/column
using `encoding` (the same encoding is applied to each column to ensure consistency).

Encoding is done internally using [`symbolize`](@ref), which typically runs a sliding
window (of width dictated by the given `encoding`) across each variable, encoding
each window to a unique integer.
"""
struct PerVariableEncoding{N} <: Discretization{N}
    outcome_spaces::NTuple{N, OutcomeSpace}
    function PerVariableEncoding(outcome_spaces::NTuple{N, OutcomeSpace}) where N
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

function PerVariableEncoding(o::OutcomeSpace)
    return PerVariableEncoding((o,))
end

function encode(encoding::PerVariableEncoding{1}, x::Vararg{<:Any, 1})
    e = first(encoding.outcome_spaces)
    x̂ = symbolize(e, first(x))
    return x̂::Vector{<:Integer}
end

function encode(encoding::PerVariableEncoding{1}, x::Tuple)
    e = first(encoding.outcome_spaces)
    x̂ = map(xᵢ -> symbolize(e, xᵢ), x)
    N = length(x)
    return x̂::NTuple{N, Vector{<:Integer}}
end

function encode(encoding::PerVariableEncoding{1}, x::AbstractStateSpaceSet)
    return encode(encoding, columns(x))
end
