using RecurrenceAnalysis: RecurrenceMatrix, JointRecurrenceMatrix
using Distances: Euclidean

export MCR
export mcr

"""
    MCR <: AssociationMeasure
    MCR(; r, metric = Euclidean())

An association measure based on mean conditional probabilities of recurrence
(MCR) introduced by Romano et al. (2007)[^Romano2007].

`r` is  mandatory keyword which specifies the recurrence threshold when constructing
recurrence matrices. It can be instance of
any subtype of `AbstractRecurrenceType` from
[RecurrenceAnalysis.jl](https://juliadynamics.github.io/RecurrenceAnalysis.jl/stable/).
To use any `r` that is not a real number, you have to do `using RecurrenceAnalysis` first.
The `metric` is any valid metric
from [Distances.jl](https://github.com/JuliaStats/Distances.jl).

## Description

For input variables `X` and `Y`, the conditional probability of recurrence
is defined as

```math
M(X | Y) = \\dfrac{1}{N} \\sum_{i=1}^N p(\\bf{y_i} | \\bf{x_i}) =
\\dfrac{1}{N} \\sum_{i=1}^N \\dfrac{\\sum_{i=1}^N J_{R_{i, j}}^{X, Y}}{\\sum_{i=1}^N R_{i, j}^X},
```

where ``R_{i, j}^X`` is the recurrence matrix and ``J_{R_{i, j}}^{X, Y}`` is the joint
recurrence matrix, constructed using the given `metric`. The measure ``M(Y | X)`` is
defined analogously.

Romano et al. (2007)'s interpretation of this quantity is that if `X` drives `Y`, then
`M(X|Y) > M(Y|X)`, if `Y` drives `X`, then `M(Y|X) > M(X|Y)`, and if coupling is symmetric,
 then `M(Y|X) = M(X|Y)`.

## Input data

`X` and `Y` can be either both univariate timeseries, or both multivariate
[`StateSpaceSet`](@ref)s.

[^Romano2007]:
    Romano, M. C., Thiel, M., Kurths, J., & Grebogi, C. (2007). Estimation of the
    direction of the coupling by conditional probabilities of recurrence. Physical Review
    E, 76(3), 036211.
"""
Base.@kwdef struct MCR{R, M} <: AssociationMeasure
    r::R
    metric::M = Euclidean()
end

function estimate(measure::MCR, x, y)
    (; r, metric) = measure
    N = length(x)
    @assert length(x) == length(y)
    jy = RecurrenceMatrix(y, r; metric)
    jxy = JointRecurrenceMatrix(x, y, r; metric)

    rp = 0.0
    for j = 1:N
        rp += @views sum(jxy[:, j]) / sum(jy[:, j])
    end
    return rp / N
end

# The
function estimate(measure::MCR, x, y, z)
    (; r, metric) = measure
    N = length(x)
    @assert length(x) == length(y)
    jy = RecurrenceMatrix(y, r; metric)
    jxy = JointRecurrenceMatrix(x, y, r; metric)
    jyz = JointRecurrenceMatrix(y, z, r; metric)
    jx_yz = JointRecurrenceMatrix(x, StateSpaceSet(y, z), r; metric)

    rp_x_y = 0.0
    rp_x_yz = 0.0
    for j = 1:N
        rp_x_y += @views sum(jxy[:, j]) / sum(jy[:, j])
        rp_x_yz += @views sum(jx_yz[:, j]) / sum(jyz[:, j])
    end
    rp_x_y /= N
    rp_x_yz /= N
    # We drop the minus sign, because we want to use negative values
    # as a criterion for rejecting dependence
    ΔMCR = -(rp_x_y - rp_x_yz)
    return ΔMCR
end

# For compatibility with causal graph and independence testing API
estimate(r::MCR, est::Nothing, x, y) = estimate(r, x, y)

"""
    mcr(m::MCR, x, y)

Compute the association between `x` and `y` based on conditional probabilities of
recurrence using the given [`MCR`](@ref) `measure`, where `x` and `y` can be either
univariate timeseries or multivariate [`StateSpaceSet`](@ref)s.
"""
mcr(args...) = estimate(args...)
