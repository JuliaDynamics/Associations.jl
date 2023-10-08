import ComplexityMeasures: codify
import ComplexityMeasures: counts
import ComplexityMeasures: probabilities

# If only one encoding is given, apply same encoding to all points
function counts(discretization::CodifyVariables{1}, x::Vararg{ArrayOrStateSpaceSet, N}) where N
    o = first(discretization.outcome_spaces)
    x̂ = (codify(o, xₖ) for xₖ in x)
    return counts(x̂...)
end

# `CodifyVariables{1}` is equivalent to `OutcomeSpace`
function probabilities(discretization::CodifyVariables{1}, x::Vararg{ArrayOrStateSpaceSet, N}) where N
    cts = counts(discretization, x...)
    return probabilities(RelativeAmount(), cts)
end
function probabilities(discretization::OutcomeSpace, x::Vararg{ArrayOrStateSpaceSet, N}) where N
    return probabilities(CodifyVariables(discretization), x...)
end
