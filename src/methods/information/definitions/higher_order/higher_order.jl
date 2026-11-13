export HigherOrderInformationMeasure

"""
    HigherOrderInformationMeasure <: MultivariateInformationMeasure

The supertype for Shannon measures of the shared information among any number of
variables, each of which is a combination of the joint entropies of subsets of
the variables.

## Implementations

- [`TotalCorrelation`](@ref)
- [`DualTotalCorrelation`](@ref)
- [`OInformation`](@ref)
- [`SInformation`](@ref)

See also [`oinformation_gradients`](@ref).

## Inputs

Each input `xᵢ` is one variable. It can be a vector (a scalar variable) or a
`StateSpaceSet` (a group of variables treated as one multivariate variable). At
least two inputs are needed.

## Compatible estimators

- [`JointProbabilities`](@ref): plug-in estimates from one joint probability mass
    function. Every marginal is taken from the same joint table, so all subset
    entropies share one discretization. The joint table is stored densely, so its
    size is the product of the number of outcomes of each variable.
- [`EntropyDecomposition`](@ref) with a [`DiscreteInfoEstimator`](@extref ComplexityMeasures.DiscreteInfoEstimator)
    of [`Shannon`](@extref ComplexityMeasures.Shannon) entropy: each subset is
    discretized and counted separately, which avoids the dense table. Any discrete
    estimator can be used, including the bias-corrected ones (`MillerMadow`, `ChaoShen`,
    `HorvitzThompson`, `Jackknife`, `Schuermann`, `GeneralizedSchuermann`), and the
    `pest` argument takes a probabilities estimator such as `BayesianRegularization`.
    Each subset entropy is estimated from the counts of that subset's *observed* joint
    outcomes; pseudo-count and shrinkage corrections therefore act on the observed
    outcomes only, not on every cell of the product of the marginal alphabets.
- [`EntropyDecomposition`](@ref) with a [`DifferentialInfoEstimator`](@extref ComplexityMeasures.DifferentialInfoEstimator)
    of [`Shannon`](@extref ComplexityMeasures.Shannon) entropy, e.g. `Kraskov`:
    each subset entropy is estimated separately. The biases of the separate
    estimates, which grow with dimension, do not cancel in general.
- [`GaussianDecomposition`](@ref): closed-form Gaussian entropies, optionally after a
    Gaussian copula transform.
"""
abstract type HigherOrderInformationMeasure <: MultivariateInformationMeasure end

min_inputs_vars(::HigherOrderInformationMeasure) = 2
max_inputs_vars(::HigherOrderInformationMeasure) = Inf

# ----------------------------------------------------------------------------------------
# Joint entropies of subsets of the input variables.
# ----------------------------------------------------------------------------------------
# For internal use. `SubsetEntropies` wraps a function `h(S)` that returns the joint
# entropy, in the unit of the measure, of the input variables with indices in `S`.
# Results are cached, because the gradients reuse the same subsets many times.
struct SubsetEntropies{F}
    h::F
    cache::Dict{Vector{Int},Float64}
end
SubsetEntropies(h) = SubsetEntropies(h, Dict{Vector{Int},Float64}())

function (H::SubsetEntropies)(S)
    key = sort!(collect(Int, S))
    isempty(key) && return 0.0
    return get!(() -> H.h(key), H.cache, key)
end

"""
    subset_entropies(est, x...) → H::SubsetEntropies

For internal use. Return a callable `H` such that `H(S)` is the estimated joint Shannon
entropy of the inputs `x[S]`, in the unit given by `est.definition.base`.
"""
function subset_entropies end

function subset_entropies(est::JointProbabilities{<:MultivariateInformationMeasure}, x...)
    p = probabilities(est.discretization, x...)
    n = length(x)
    measure = Shannon(; base = est.definition.base)
    h(S) = information(measure, length(S) == n ? p : marginal(p; dims = S))
    return SubsetEntropies(h)
end

function subset_entropies(
        est::EntropyDecomposition{<:MultivariateInformationMeasure, <:DifferentialInfoEstimator{<:Shannon}},
        x...
    )
    X = StateSpaceSet.(x)
    modified_est = estimator_with_overridden_parameters(est.definition, est.est)
    h(S) = information(modified_est, StateSpaceSet(X[S]...))
    return SubsetEntropies(h)
end

function subset_entropies(
        est::EntropyDecomposition{<:MultivariateInformationMeasure, <:DiscreteInfoEstimator{<:Shannon}},
        x...
    )
    # One encoding per input; a multi-column input gives a multi-column encoding.
    E = StateSpaceSet.(codified_marginals(est.discretization, x...))
    modified_est = estimator_with_overridden_parameters(est.definition, est.est)
    o = UniqueElements()
    h(S) = information(modified_est, est.pest, o, StateSpaceSet(E[S]...))
    return SubsetEntropies(h)
end

# ----------------------------------------------------------------------------------------
# Each measure is a combination of subset entropies. `entropy_combination(def, H, V)`
# evaluates the measure on the variables with indices `V`.
# ----------------------------------------------------------------------------------------
function entropy_combination end

const HIGHER_ORDER_ESTIMATORS = Union{
    JointProbabilities{<:HigherOrderInformationMeasure},
    EntropyDecomposition{<:HigherOrderInformationMeasure},
}

function association(est::HIGHER_ORDER_ESTIMATORS, x...)
    return higher_order_association(est, x...)
end

function higher_order_association(est, x...)
    verify_number_of_inputs_vars(est.definition, length(x))
    H = subset_entropies(est, x...)
    return entropy_combination(est.definition, H, collect(1:length(x)))
end

include("TotalCorrelation.jl")
include("DualTotalCorrelation.jl")
include("OInformation.jl")
include("SInformation.jl")
include("oinformation_gradients.jl")
