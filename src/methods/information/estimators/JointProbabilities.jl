using ComplexityMeasures: OutcomeSpace
export JointProbabilities

"""
    JointProbabilities <: InformationMeasureEstimator
    JointProbabilities(
        definition::MultivariateInformationMeasure,
        discretization::Discretization
    )

`JointProbabilities` is a generic estimator for multivariate discrete information measures.
    
## Usage

- Use with [`association`](@ref) to compute an information measure from input data.

## Description

It first encodes the input data according to the given `discretization`, then constructs 
`probs`, a multidimensional [`Probabilities`](@ref) instance. Finally, `probs` are 
forwarded to a [`PlugIn`](@ref) estimator, which computes the measure according to 
`definition`.

# Compatible encoding schemes

-  [`CodifyVariables`](@ref) (encode each *variable*/column of the input data independently by 
    applying an encoding in a sliding window over each input variable).  
- [`CodifyPoints`](@ref) (encode each *point*/column of the input data)

Works for any [`OutcomeSpace`](@ref) that implements [`codify`](@ref).

!!! note "Joint probabilities vs decomposition methods"

    Using [`JointProbabilities`](@ref) to compute an information measure, 
    e.g. conditional mutual estimation,
    is typically slower than other dedicated estimation procedures like [`EntropyDecomposition`](@ref).
    The reason is that measures such as [`CMIShannon`](@ref) can be formulated as a
    sum of four entropies, which can be estimated individually and summed afterwards. 
    This decomposition is fast because because we avoid *explicitly* estimating the entire joint pmf, 
    which demands many extra calculation steps, However, the decomposition is biased, 
    because it fails to fully take into consideration the joint relationships between the variables.
    Pick your estimator according to your needs.

See also: [`Counts`](@ref), [`Probabilities`](@ref), [`ProbabilitiesEstimator`](@ref),
[`OutcomeSpace`](@ref), [`DiscreteInfoEstimator`](@ref).
"""
struct JointProbabilities{M <: MultivariateInformationMeasure, O, P} <: MultivariateInformationMeasureEstimator{M}
    definition::M # API from complexity measures: definition must be the first field of the infoestimator.
    discretization::O
    pest::P # Not exposed to user for now.

    function JointProbabilities(def::M, disc::D, pest = RelativeAmount()) where {M, D}
        new{M, D, typeof(pest)}(def, disc, pest)
    end
end

function association(est::JointProbabilities{MultivariateInformationMeasure}, x...)
    throw(ArgumentError("`JointProbabilities` not implemented for `$(typeof(est.definition).name.name)`"))
end