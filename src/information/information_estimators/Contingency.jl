export Contingency

"""
    Contingency <: DiscreteInfoEstimator
    Contingency(definition::MultivariateInformationMeasure,
        pest::ProbabilitiesEstimator, o::OutcomeSpace)

`Contingency` is a generic estimator for discrete information measures that first
discretizes/encodes the input data according to the given [`OutcomeSpace`](@ref) `o`, then
constructs a contingency table of the required dimensionality (a [`Counts`](@ref) instance),
then constructs a multidimensional probability mass function (a [`Probabilities`](@ref)
instance) from those counts using the provided probabilities estimator `pest`.

Works for any outcome space that implements [`symbolize`](@ref).

See also: [`Counts`](@ref), [`Probabilities]](@ref), [`ProbabilitiesEstimator`](@ref),
[`OutcomeSpace`](@ref), [`DiscreteInfoEstimator`](@ref).
"""
Base.@kwdef struct Contingency{M <: MultivariateInformationMeasure, O, P} <: DiscreteInfoEstimator
    definition::M # API from complexity measures: definition must be first field.
    o::O
    pest::P
end
