using ComplexityMeasures: OutcomeSpace
export JointProbabilities

"""
    JointProbabilities <: InformationMeasureEstimator
    JointProbabilities(
        definition::MultivariateInformationMeasure,
        discretization::Discretization
    )

`JointProbabilities` is a generic estimator for discrete information measures that first
discretizes/encodes the input data according to the `discretization` (either a 
[`CodifyVariables`](@ref) in the case of time series or the more generic
[`CodifyPoints`](@ref) in other cases), then
constructs a contingency table of the required dimensionality (a [`Counts`](@ref) instance),
then constructs a multidimensional probability mass function (a [`Probabilities`](@ref)
instance) using plug-in estimation of probabilities (relative frequencies of counts).

Works for any outcome space that implements [`codify`](@ref).

See also: [`Counts`](@ref), [`Probabilities`](@ref), [`ProbabilitiesEstimator`](@ref),
[`OutcomeSpace`](@ref), [`DiscreteInfoEstimator`](@ref).
"""
struct JointProbabilities{M <: MultivariateInformationMeasure, O, P} <: MultivariateInformationMeasureEstimator{M}
    definition::M # API from complexity measures: definition must be the first field of the infoestimator.
    discretization::O
    pest::P # Not exposed for user for now.

    function JointProbabilities(def::M, disc::D, pest = RelativeAmount()) where {M, D}
        new{M, D, typeof(pest)}(def, disc, pest)
    end
end

# One method to compute them all (through joint probabilities)
const DiscretizationOrOutcomeSpace = Union{Discretization, OutcomeSpace}
function information(est::JointProbabilities{<:M, <:DiscretizationOrOutcomeSpace}, x...) where M
    probs = probabilities(est.discretization, x...)
    return information(est.definition, probs)
end
