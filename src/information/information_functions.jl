"""
    information(definition::MultivariateInformationMeasure, p::Probabilities{T, N})

Estimate an `N`-variate multivariate information measure given by `definition` 
directly from a  pre-computed joint probability mass function `p`.
"""
function information(::MultivariateInformationMeasure, p::Probabilities) end

"""
    information(est::MultivariateInformationMeasureEstimator, x, y, [z, w, ...])

Using the provided estimator `est`, estimate the information measure given by
`est.definition` from input data `x, y, ...`, where the number of inputs depend on `est`.

The docstring for [`MultivariateInformationMeasure`](@ref) lists all possible measures,
and the docstring for [`MultivariateInformationMeasureEstimator`](@ref) lists possible
estimators. You may also want to check out the
[convenience wrappers](@ref convenience_info) in the online documentation.


## Examples

```julia
using CausalityTools
using Random; rng = MersenneTwister(12345)
x, y, z = rand(rng, 100), rand(rng, 100), rand(rng, 100)


# Mutual information
# ------------------
est = EntropyDecomposition(MIShannon(), PlugIn(Shannon()), OrdinalPatterns(m=3))
information(est, x, y)

est = JointProbabilities(MIShannon(), CodifyVariables(OrdinalPatterns(m=3)))
information(est, x, y)

# Conditional mutual information
information(JointProbabilities(CMIShannon(), ValueBinning(3)), x, y, z)
```
"""
function information(::MultivariateInformationMeasureEstimator, x...) end


"""
    information(definition::BivariateInformationMeasure, d::Discretization, x, y)

Estimate the bivariate information measure according to its `definition`, from
input data `x` and `y` (where `length(x) == length(y)`).

This is done by first discretizing `x` and `y` according to the discretization scheme `d`,
then estimating a 2D probability mass function (pmf) over the distinct outcomes dictated by
`d`, and finally applying the `definition` to those probabilities.

See also: [`Discretization`](@ref), [`BivariateInformationMeasure`](@ref).
"""
function information(measure::BivariateInformationMeasure,
            pest::ProbabilitiesEstimator,
            d::Discretization,
            x, y)

end
function information(measure::BivariateInformationMeasure, o::OutcomeSpace, x...)
    return information(measure, RelativeAmount(), o, x...)
end

# A generic implementation suffices. It dispatches to measure-specific functions
# that only depend on two or more sets of `Probabilities`.
function information(measure::MultivariateInformationMeasure, est::ProbabilitiesEstimator, o::OutcomeSpace, x...)
    cts_xs = map(xᵢ -> first(allcounts_and_outcomes(o, xᵢ)), x)
    pmfs = map(cts -> probabilities(est, cts), cts_xs)
    return information(measure, pmfs...)
end

# For the divergences/distances, is it faster to compute the joint and derive the marginal from the joint, 
# or use allcounts to compute the marginals directly? We stick with the joint for now...
# function information(measure::BivariateInformationMeasure, est::ProbabilitiesEstimator, o::OutcomeSpace, x, y)
#     cts_x, outs_x = allcounts_and_outcomes(o, x)
#     cts_y, outs_y = allcounts_and_outcomes(o, y)
#     px = probabilities(est, cts_x)
#     py = probabilities(est, cts_y)
#     return information(measure, px, py)
# end


function size_match(measure::BivariateInformationMeasure, px::Probabilities, py::Probabilities)
    size(px) == size(py) || throw(DimensionMismatch("px and py must have the same size"))
end