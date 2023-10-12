"""
    information(definition::MultivariateInformationMeasure, p::Probabilities{T, N})

Estimate an `N`-variate multivariate information measure given by `definition` 
directly from a  pre-computed joint probability mass function `p`.
"""
function information(::MultivariateInformationMeasure, p::Probabilities) end

"""
    information(est::MultivariateInformationMeasureEstimator, x...)

Estimate the information measure given by `est.definition` from the input data `x`,
where `length(x) ≥ 2`.

## Estimators

- [`JointProbabilities`](@ref). 
- [`EntropyDecomposition`](@ref). 
- [`MIDecomposition`](@ref).

When `length(x) == 2`, any of the following estimators can be used:

- Any [`MutualInformationEstimator`](@ref)

When `length(x) == 3`, any of the following estimators can be used:

- Any [`ConditionalMutualInformationEstimator`](@ref)


## Examples

Here's three different ways 
```julia
using CausalityTools
est = JointProbabilities(CMIShannon(), ValueBinning(3))
x, y, z = rand(100), rand(100), rand(100)
information(est, x, y, z)
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
    cts_xs = map(xᵢ -> allcounts(o, xᵢ), x)
    pmfs = map(cts -> probabilities(est, cts), cts_xs)
    return information(measure, pmfs...)
end

# For the divergences/distances, is it faster to compute the joint and derive the marginal from the joint, 
# or use allcounts to compute the marginals directly? We stick with the joint for now...
# function information(measure::BivariateInformationMeasure, est::ProbabilitiesEstimator, o::OutcomeSpace, x, y)
#     cts_x = allcounts(o, x)
#     cts_y = allcounts(o, y)
#     px = probabilities(est, cts_x)
#     py = probabilities(est, cts_y)
#     return information(measure, px, py)
# end


function size_match(measure::BivariateInformationMeasure, px::Probabilities, py::Probabilities)
    size(px) == size(py) || throw(DimensionMismatch("px and py must have the same size"))
end
