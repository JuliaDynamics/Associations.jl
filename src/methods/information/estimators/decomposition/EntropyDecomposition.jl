
export EntropyDecomposition

"""
    EntropyDecomposition(definition::MultivariateInformationMeasure, 
        est::DifferentialInfoEstimator)
    EntropyDecomposition(definition::MultivariateInformationMeasure,
        est::DiscreteInfoEstimator,
        discretization::CodifyVariables{<:OutcomeSpace},
        pest::ProbabilitiesEstimator = RelativeAmount())

Estimate the multivariate information measure specified by `definition` by rewriting
its formula into some combination of entropy terms. 

## Estimation 

`EntropyDecomposition` allows using any 
[`InformationMeasureEstimators`s](@extref ComplexityMeasures.InformationMeasureEstimator) from 
[ComplexityMeasures.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/complexitymeasures/stable/)
to estimate multivariate information measures.

Computations are done by computing individual entropy terms using `est`, then combining them according to 
`definition` to form the final estimate. 

## Bias 

Estimating the `definition` by decomposition into a combination of entropy terms,
which are estimated independently, will in general be more biased than when using a
dedicated estimator. One reason is that this decomposition may miss out on crucial
information in the joint space. To remedy this, dedicated information measure 
estimators typically derive the marginal estimates by first considering the joint
space, and then does some clever trick to eliminate the bias that is introduced
through a naive decomposition. Unless specified below, no bias correction is 
applied for `EntropyDecomposition`.

## Usage

- Use with [`association`](@ref) to compute a [`MultivariateInformationMeasure`](@ref)
    from input data: [`association`](@ref)`(est::EntropyDecomposition, x...)`.
- Use with some [`IndependenceTest`](@ref) to test for independence between variables.

## Description

If `est` is a [`DifferentialInfoEstimator`](@extref ComplexityMeasures.DifferentialInfoEstimator), then `discretization` and `pest` 
are ignored. 
## Differential estimation

If using the first signature, any compatible [`DifferentialInfoEstimator`](@extref ComplexityMeasures.DifferentialInfoEstimator) can be 
used.
[`MITsallisMartin`](@ref) can be estimated using a decomposition into entropy 
terms using [`EntropyDecomposition`](@ref). This is done by using estimators from 
[ComplexityMeasures.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/complexitymeasures/stable/). We can use any compatible 
[`InformationMeasureEstimator`](@extref ComplexityMeasures.InformationMeasureEstimator)
that can estimate differential [`Tsallis`](@extref ComplexityMeasures.Tsallis) entropy from
[ComplexityMeasures.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/complexitymeasures/stable/).

## Discrete estimation

If `est` is a [`DiscreteInfoEstimator`](@extref ComplexityMeasures.DiscreteInfoEstimator), then 
`discretization` and a probabilities estimator `pest` must also be provided (default to `RelativeAmount`, 
which uses naive plug-in probabilities). This will always discretize the input data 
per variable/column, and the encode the discretized column into integers using [`codify`](@ref).

The following [`OutcomeSpace`s](@extref ComplexityMeasures.OutcomeSpace) can be used for discretisation. 
Note that not all outcome spaces will work with all measures.

| Estimator                                                                       | Principle                     |
| :------------------------------------------------------------------------------ | :---------------------------- |
| [`UniqueElements`](@extref ComplexityMeasures.UniqueElements)                   | Count of unique elements      |
| [`ValueBinning`](@extref ComplexityMeasures.ValueBinning)                       | Binning (histogram)           |
| [`OrdinalPatterns`](@extref ComplexityMeasures.OrdinalPatterns)                 | Ordinal patterns              |
| [`Dispersion`](@extref ComplexityMeasures.Dispersion)                           | Dispersion patterns           |
| [`BubbleSortSwaps`](@extref ComplexityMeasures.BubbleSortSwaps)                 | Sorting complexity            |
| [`CosineSimilarityBinning`](@extref ComplexityMeasures.CosineSimilarityBinning) | Cosine similarities histogram |


## Handling of overlapping parameters

If there are overlapping parameters between the measure to be estimated, and the
lower-level decomposed measures, then the top-level measure parameter takes precedence.
For example, if we want to estimate `CMIShannon(base = 2)` through a decomposition 
of entropies using the `Kraskov(Shannon(base = ℯ))` Shannon entropy estimator, then
`base = 2` is used.

!!! info 
    Not all measures have the property that they can be decomposed into more fundamental
    information theoretic quantities. For example, [`MITsallisMartin`](@ref) *can* be 
    decomposed into a combination of marginal entropies, while [`MIRenyiSarbu`](@ref)
    cannot. An error will be thrown if decomposition is not possible.

## Discrete entropy decomposition 

The second signature is for discrete estimation using [`DiscreteInfoEstimator`](@extref ComplexityMeasures.DiscreteInfoEstimator)s,
for example [`PlugIn`](@extref ComplexityMeasures.PlugIn). The given `discretization` scheme (typically an 
[`OutcomeSpace`](@extref ComplexityMeasures.OutcomeSpace)) controls how the joint/marginals are discretized, and the
probabilities estimator `pest` controls how probabilities are estimated from counts.

!!! note "Bias"
    Like for [`DifferentialInfoEstimator`](@extref ComplexityMeasures.DifferentialInfoEstimator), using a dedicated estimator 
    for the measure in question will be more reliable than using a decomposition
    estimate. Here's how different `discretization`s are applied:

    - [`ValueBinning`](@extref ComplexityMeasures.ValueBinning). Bin visitation frequencies are counted in the joint space
        `XY`, then marginal visitations are obtained from the joint bin visits.
        This behaviour is the same for both [`FixedRectangularBinning`](@extref ComplexityMeasures.FixedRectangularBinning) and
        [`RectangularBinning`](@extref ComplexityMeasures.RectangularBinning) (which adapts the grid to the data).
        When using [`FixedRectangularBinning`](@extref ComplexityMeasures.FixedRectangularBinning), the range along the first dimension
        is used as a template for all other dimensions. This is a bit slower than naively 
        binning each marginal, but lessens bias.
    - [`OrdinalPatterns`](@extref ComplexityMeasures.OrdinalPatterns). Each timeseries is separately [`codify`](@ref)-ed
        according to its ordinal pattern (no bias correction).
    - [`Dispersion`](@extref ComplexityMeasures.Dispersion). Each timeseries is separately [`codify`](@ref)-ed according
        to its dispersion pattern  (no bias correction).

## Examples

- [Example 1](@ref example_MIShannon_EntropyDecomposition_Jackknife_ValueBinning):
    [`MIShannon`](@ref) estimation using decomposition into discrete [`Shannon`](@extref ComplexityMeasures.Shannon)
    entropy estimated using [`CodifyVariables`](@ref) with [`ValueBinning`](@extref ComplexityMeasures.ValueBinning).
- [Example 2](@ref example_MIShannon_EntropyDecomposition_BubbleSortSwaps):
    [`MIShannon`](@ref) estimation using decomposition into discrete [`Shannon`](@extref ComplexityMeasures.Shannon)
    entropy estimated using [`CodifyVariables`](@ref) with [`BubbleSortSwaps`](@extref ComplexityMeasures.BubbleSortSwaps).
- [Example 3](@ref example_MIShannon_EntropyDecomposition_Kraskov):
    [`MIShannon`](@ref) estimation using decomposition into differental [`Shannon`](@extref ComplexityMeasures.Shannon)
    entropy estimated using the [`Kraskov`](@extref ComplexityMeasures.Kraskov) estimator.

See also: [`MutualInformationEstimator`](@ref), [`MultivariateInformationMeasure`](@ref).
"""
struct EntropyDecomposition{
    M<:MultivariateInformationMeasure,
    E<:InformationMeasureEstimator,
    D<:Union{Discretization,Nothing},
    P<:Union{ProbabilitiesEstimator,Nothing}
} <: DecompositionEstimator{M}
    definition::M # extend API from complexity measures: definition must be the first field of the info estimator.
    est::E # The estimator + measure which `definition` is decomposed into.
    discretization::D # `Nothing` if `est` is a `DifferentialInfoEstimator`.
    pest::P # `Nothing` if `est` is a `DifferentialInfoEstimator`.


    function EntropyDecomposition(
        definition::MultivariateInformationMeasure,
        est::DifferentialInfoEstimator)
        M = typeof(definition)
        E = typeof(est)
        verify_decomposition_entropy_type(definition, est)
        return new{M,E,Nothing,Nothing}(definition, est, nothing, nothing)
    end

    function EntropyDecomposition(
        definition::MultivariateInformationMeasure,
        est::DiscreteInfoEstimator,
        discretization::D,
        pest::ProbabilitiesEstimator=RelativeAmount(),
    ) where {D}
        M = typeof(definition)
        E = typeof(est)
        P = typeof(pest)
        verify_decomposition_entropy_type(definition, est)

        return new{M,E,D,P}(definition, est, discretization, pest)
    end
end

# For internal use.
"""
    verify_decomposition_entropy_type(
        definition::MultivariateInformationMeasure, 
        est::Union{DiscreteInfoEstimator, DifferentialInfoEstimator}
    )

Check that we can actually decompose the `definition` into `est.definition`. The 
default is to do nothing. Certain definitions  may override (e.g. `CMIRenyiJizba` does so).
"""
function verify_decomposition_entropy_type(
    definition::MultivariateInformationMeasure,
    est::Union{DiscreteInfoEstimator,DifferentialInfoEstimator})
end


# ----------------------------------------------------------------------------------------
# Custom pretty printing for discrete entropy estimators, since it has more field.
# ----------------------------------------------------------------------------------------
function summary_strings(est::EntropyDecomposition{<:M,<:DiscreteInfoEstimator}) where M
    return [
        "Measure to be decomposed",
        "Estimator for decomposed components",
        "Discretization",
        "Probabilities estimator"
    ]
end

function summary_types(est::EntropyDecomposition{<:M,<:DiscreteInfoEstimator}) where M
    return [
        typeof(est.definition),
        typeof(est.est),
        typeof(est.discretization),
        typeof(est.pest)
    ]
end

function measure_colors(est::EntropyDecomposition{<:M,<:DiscreteInfoEstimator}) where M
    return [
        :light_red,
        :light_green,
        :light_blue,
        :light_yellow,
    ]
end

function info_colors(est::EntropyDecomposition{<:M,<:DiscreteInfoEstimator}) where M
    return [
        :red,
        :green,
        :blue,
        :yellow,
    ]
end
