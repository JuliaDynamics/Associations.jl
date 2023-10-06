import ComplexityMeasures: probabilities
export marginal

# ##########################################################################################
# Probabilities API.
# The following code extends the functionality of ComplexityMeasures.jl for multiple
# input variables (ComplexityMeasures.jl only deals with single-variable estimation)
# ##########################################################################################
function probabilities(est::RelativeAmount, c::Counts{<:Integer, N}) where N
    return Probabilities(c)
end

function probabilities(d::Discretization, x::Vararg{Any, N}) where N
    cts = counts(d, x...)
    return probabilities(RelativeAmount(), cts)
end

# Not providing any discretization defaults to `RelativeAmount` estimation.
function probabilities(x::Vararg{VectorOrStateSpaceSet, N}) where N
    cts = counts(UniqueElements(), x...)
    return probabilities(RelativeAmount(), cts)
end

# TODO: preserve axis labels
"""
    marginal(p::Probabilities; dims = 1:ndims(p))
    marginal(c::Counts; dims = 1:ndims(p))

Given a set of counts `c` (a contingency table), or a multivariate probability mass
function `p`, return the marginal counts/probabilities along the given `dims`.
"""
function marginal(p::Probabilities; dims = 1:ndims(p))
    alldims = 1:ndims(p)
    reduce_dims = (setdiff(alldims, dims)...,)
    marginal = dropdims(sum(p.p, dims = reduce_dims), dims = reduce_dims)
    return Probabilities(marginal)
end
