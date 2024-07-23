"""
    information(definition::MultivariateInformationMeasure, p::Probabilities{T, N})

Estimate an `N`-variate multivariate information measure given by `definition` 
directly from a  pre-computed joint probability mass function `p`.
"""
function information(::MultivariateInformationMeasure, p::Probabilities) end

function size_match(measure::BivariateInformationMeasure, px::Probabilities, py::Probabilities)
    size(px) == size(py) || throw(DimensionMismatch("px and py must have the same size"))
end
