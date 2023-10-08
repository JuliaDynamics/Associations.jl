using ComplexityMeasures
import ComplexityMeasures: codify
import ComplexityMeasures: OutcomeSpace

using DelayEmbeddings: embed
export CodifyVariables
export codify

# TODO: implement this Generically for `Encodings` too (will require type-parameterized
# number of elements for the Encodings).

"""
    CodifyVariables <: Discretization
    CodifyVariables(outcome_spaces::OutcomeSpace)

Given multiple dataset `xs::StateSpaceSet`, [`encode`](@ref) the `i`-th variable/column
using `encoding` (the same encoding is applied to each column to ensure consistency).

Encoding is done internally using [`codify`](@ref), which typically runs a sliding
window (of width dictated by the given `encoding`) across each variable, encoding
each window to a unique integer.
"""
struct CodifyVariables{N} <: Discretization{N}
    outcome_spaces::NTuple{N, OutcomeSpace}
    function CodifyVariables(outcome_spaces::NTuple{N, OutcomeSpace}) where N
        if N > 1
            s = "It is currently only possible to use the same `OutcomeSpace` for all " *
                "variables. Got $N different encodings"
            throw(ArgumentError(s))
        end
        new{N}(outcome_spaces)
    end
end

function CodifyVariables(o::OutcomeSpace)
    return CodifyVariables((o,))
end

function codify(encoding::CodifyVariables{1}, x::Vararg{Any, 1})
    e = first(encoding.outcome_spaces)
    x̂ = ComplexityMeasures.codify(e, first(x))
    return x̂::Vector{<:Integer}
end

function codify(encoding::CodifyVariables{1}, x::Tuple)
    e = first(encoding.outcome_spaces)
    x̂ = map(xᵢ -> ComplexityMeasures.codify(e, xᵢ), x)
    N = length(x)
    return x̂::NTuple{N, Vector{<:Integer}}
end

function codify(encoding::CodifyVariables{1}, x::AbstractStateSpaceSet)
    return codify(encoding, columns(x))
end
