
export EntropyDecomposition

"""
    EntropyDecomposition(definition::MultivariateInformationMeasure, 
        est::DifferentialInfoEstimator)
    EntropyDecomposition(definition::MultivariateInformationMeasure,
        est::DiscreteInfoEstimator,
        discretization::OutcomeSpace,
        pest::ProbabilitiesEstimator = RelativeAmount())

If `est` is a [`DifferentialInfoEstimator`](@ref), then `discretization` and `pest` 
are ignored. If `est` is a [`DiscreteInfoEstimator`](@ref), then `discretization` and a
probabilities estimator `pest` must also be provided (default to `RelativeAmount`, 
which uses naive plug-in probabilities).

## Usage

- [`information`](@ref)`(est::EntropyDecomposition, x...)`.

See also: [`MutualInformationEstimator`](@ref), [`MultivariateInformationMeasure`](@ref).
"""
struct EntropyDecomposition{M <: MultivariateInformationMeasure, E <: InformationMeasureEstimator, D, P} <: InformationMeasureEstimator{M}
    definition::M # extend API from complexity measures: definition must be the first field of the info estimator.
    est::E # The estimator + measure which `definition` is decomposed into.
    discretization::D # `Nothing` if `est` is a `DifferentialInfoEstimator`.
    pest::P # `Nothing` if `est` is a `DifferentialInfoEstimator`.


    function EntropyDecomposition(
        definition::MultivariateInformationMeasure, 
        est::DifferentialInfoEstimator)
        M = typeof(definition)
        E = typeof(est)
        return new{M, E, Nothing, Nothing}(definition, est, nothing, nothing)
    end

    function EntropyDecomposition(
            definition::MultivariateInformationMeasure, 
            est::DiscreteInfoEstimator, 
            discretization::D,
            pest::ProbabilitiesEstimator = RelativeAmount(),
        ) where {D}
        M = typeof(definition)
        E = typeof(est)
        P = typeof(pest)
    
        return new{M, E, D, P}(definition, est, discretization, pest)
    end

end


