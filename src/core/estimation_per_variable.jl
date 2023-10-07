import ComplexityMeasures: codify
import ComplexityMeasures: counts
import ComplexityMeasures: probabilities

# If only one encoding is given, apply same encoding to all points
function counts(discretization::CodifyVariables{1}, x::Vararg{ArrayOrStateSpaceSet, N}) where N
    o = first(discretization.outcome_spaces)
    x̂ = (codify(o, xₖ) for xₖ in x)
    return counts(x̂...)
end

# # If a single outcome space is given, then codify each variable according to the same outcome space.
# function counts(discretization::OutcomeSpace, x::Vararg{ArrayOrStateSpaceSet, 1}) where N
#     return ComplexityMeasures.counts(discretization, x)
# end


# # If a single outcome space is given, then codify each variable according to the same outcome space.
# function counts(discretization::OutcomeSpace, x::Vararg{ArrayOrStateSpaceSet, N}) where N
#     return counts(CodifyVariables(discretization), x)
# end

# function probabilities(discretization::CodifyVariables{1}, x::Vararg{ArrayOrStateSpaceSet, N}) where N
#     cts = counts(discretization, x...)
#     return probabilities(RelativeAmount(), cts)
# end
# # Equivalent to the above.
# function probabilities(o::OutcomeSpace, x::Vararg{Any, N}) where N
#     cts = counts(o, x...)
#     return probabilities(RelativeAmount(), cts)
# end
