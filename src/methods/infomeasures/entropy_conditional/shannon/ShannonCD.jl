export ShannonCD



"""
    ShannonCD <: ConditionalDifferentialEntropyEstimator
    ShannonCD(est::DifferentialEntropyEstimator)
    ShannonCD(est::ProbabilitiesEstimator)

`ShannonCD` is a generic estimator that computes the discrete or differential conditional
Shannon entropy, depending on `est`.

## Description

### Discrete Shannon conditional entropy

If `est` is a [`ConditionalDifferentialEntropyEstimator`](@ref), then compute the discrete conditional
entropy by first approximating probabilities, and plugging them into the formula below.
The `ShannonCD` estimator thus can use any [`ProbabilitiesEstimator`](@ref) to compute a
plug-in estimate of``H(X | Y) = H(X,Y) - H(X)``.

!!! warn "Common outcome space"
    Input variables must share [`outcome_space`](@ref) for probability-based estimators to
    yield meaningful estimates. For example, if using [`ValueHistogram`](@ref), make
    sure to use a [`FixedRectangularBinning`](@ref).

## Continuous/differential Shannon conditional entropy

If `est` is an [`DifferentialEntropyEstimator`](@ref) then the differential conditional entropy
``h(X | Y) = -\\int h(X,Y) - h(X)``, where ``h(\\cdot)`` is the differential entropy.

"""
struct ShannonCD{E} <: ConditionalDifferentialEntropyEstimator
    est::E
end

function entropy_conditional(e::Renyi, est::ShannonCD, x, y)
    e.q == 1 || throw(ArgumentError(
        "Renyi entropy with q = $(e.q) not implemented for $(typeof(est)) estimators"
    ))
    entropy(e, est.est, Dataset(x, y)) - entropy(e, est.est, Dataset(x))
end
