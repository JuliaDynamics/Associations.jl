# ##########################################################################################
# Probabilities API.
# The following code extends the functionality of ComplexityMeasures.jl for multiple
# input data (ComplexityMeasures.jl only deals with single-variable estimation)
# ##########################################################################################
function probabilities(est::RelativeAmount, c::Counts{<:Integer, N}) where N
    return Probabilities(c)
end

function probabilities(d::Discretization, x::Vararg{Any, N}) where N
    cts = counts(d, x...)
    return probabilities(RelativeAmount(), cts)
end
