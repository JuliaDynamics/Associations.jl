
export EntropyDecomposition

"""
    EntropyDecomposition(definition::MultivariateInformationMeasure, 
        est::DifferentialInfoEstimator)
    EntropyDecomposition(definition::MultivariateInformationMeasure,
        est::DiscreteInfoEstimator,
        discretization::OutcomeSpace,
        pest::ProbabilitiesEstimator = RelativeAmount())

Estimate the multivariate information measure specified by `definition` by rewriting
its formula into some combination of entropy terms. 

The entropy terms are estimated using `est`, and then combined to form the final 
estimate of `definition`. No bias correction is applied.
If `est` is a [`DifferentialInfoEstimator`](@ref), then `discretization` and `pest` 
are ignored. If `est` is a [`DiscreteInfoEstimator`](@ref), then `discretization` and a
probabilities estimator `pest` must also be provided (default to `RelativeAmount`, 
which uses naive plug-in probabilities).

## Bias 

Estimating the `definition` by decomposition into a combination of entropy terms,
which are estimated independently, will in general be more biased than when using a
dedicated estimator. One reason is that this decomposition may miss out on crucial
information in the joint space. To remedy this, dedicated information measure 
estimators typically derive the marginal estimates by first considering the joint
space, and then does some clever trick to eliminate the bias that is introduced
through a naive decomposition. Unless specified below, no bias correction is 
applied for `EntropyDecomposition`.


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

The second signature is for discrete estimation using [`DiscreteInfoEstimator`](@ref)s,
for example [`PlugIn`](@ref). The given `discretization` scheme (typically an 
[`OutcomeSpace`](@ref)) controls how the joint/marginals are discretized, and the
probabilities estimator `pest` controls how probabilities are estimated from counts.

!!! note "Bias"
    Like for [`DifferentialDecomposition`](@ref), using a dedicated estimator 
    for the measure in question will be more reliable than using a decomposition
    estimate. Here's how different `discretization`s are applied:

    - [`ValueBinning`](@ref). Bin visitation frequencies are counted in the joint space
        `XY`, then marginal visitations are obtained from the joint bin visits.
        This behaviour is the same for both [`FixedRectangularBinning`](@ref) and
        [`RectangularBinning`](@ref) (which adapts the grid to the data).
        When using [`FixedRectangularBinning`](@ref), the range along the first dimension
        is used as a template for all other dimensions. This is a bit slower than naively 
        binning each marginal, but lessens bias.
    - [`OrdinalPatterns`](@ref). Each timeseries is separately [`codify`](@ref)-ed
        according to its ordinal pattern (no bias correction).
    - [`Dispersion`](@ref). Each timeseries is separately [`codify`](@ref)-ed according
        to its dispersion pattern  (no bias correction).

## Usage

- [`information`](@ref)`(est::EntropyDecomposition, x...)`.

See also: [`MutualInformationEstimator`](@ref), [`MultivariateInformationMeasure`](@ref).

## Examples

Both Shannon-type mutual information and conditional mutual information can be 
written as a sum of marginal entropy terms. First a discrete example for mutual 
information:

```julia
using CausalityTools
using Random; rng = MersenneTwister(1234)

x = StateSpaceSet(rand(rng, 1000000, 2))
y = StateSpaceSet(rand(rng, 1000000, 2))
# Compute Shannon mutual information by discretizing each marginal column-wise
# (per variable) using length-`3` ordinal patterns.
est = EntropyDecomposition(MIShannon(), PlugIn(Shannon()), OrdinalPatterns(m=3))
information(est, x, y) # should be close to 0
```

Here, we estimate Shannon-type conditional mutual information using the `ZhuSingh`
entropy estimator.

```julia
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 100000)
y = rand(rng, 100000) .+ x
z = rand(rng, 100000) .+ y

est = EntropyDecomposition(CMIShannon(), ZhuSingh(k = 3))
information(est, x, z, y) # should be near 0 (and can be negative)
```
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

# ----------------------------------------------------------------------------------------
# Pretty printing
# ----------------------------------------------------------------------------------------
function summary_strings(est::EntropyDecomposition{<:M, <:DifferentialInfoEstimator}) where M
    return [
        "Measure to be decomposed",
        "Estimator for decomposed components",
    ]
end

function summary_strings(est::EntropyDecomposition{<:M, <:DiscreteInfoEstimator}) where M
    return [
        "Measure to be decomposed",
        "Estimator for decomposed components",
        "Discretization",
        "Probabilities estimator"
    ]
end

function summary_types(est::EntropyDecomposition{<:M, <:DifferentialInfoEstimator}) where M
    return [
        typeof(est.definition),
        typeof(est.est),
    ]
end

function summary_types(est::EntropyDecomposition{<:M, <:DiscreteInfoEstimator}) where M
    return [
        typeof(est.definition),
        typeof(est.est),
        typeof(est.discretization),
        typeof(est.pest)
    ]
end

function measure_colors(est::EntropyDecomposition{<:M, <:DifferentialInfoEstimator}) where M
    return [
        :light_red,
        :light_green,
    ]
end

function measure_colors(est::EntropyDecomposition{<:M, <:DiscreteInfoEstimator}) where M
    return [
        :light_red,
        :light_green,
        :light_blue,
        :light_yellow,
    ]
end

function info_colors(est::EntropyDecomposition{<:M, <:DifferentialInfoEstimator}) where M
    return [
        :red,
        :green,
    ]
end

function info_colors(est::EntropyDecomposition{<:M, <:DiscreteInfoEstimator}) where M
    return [
        :red,
        :green,
        :blue,
        :yellow,
    ]
end

function Base.show(io::IO, est::EntropyDecomposition)
    types = summary_types(est)
    strs = summary_strings(est)
    measurecolors = measure_colors(est)
    infocolors = info_colors(est)
    n = maximum(length.(strs))

    spaces_needed = [n - length(s) for s in strs] 
    spaced_strs = [strs[i] * repeat(" ", spaces_needed[i]) for i in eachindex(strs)]
    ctx = IOContext(io, :color => true)
    printstyled(ctx,  "EntropyDecomposition estimator\n\n", color=:bold)
    d = decomposition_string(est.definition, est)
        printstyled(ctx,  "  Formula: $(d)\n\n", color=:light_grey)
    indent = " "
    for i in eachindex(strs)
        printstyled(ctx, "$(indent)$(spaced_strs[i]): ", color=infocolors[i])
        printstyled(ctx, string(types[i]), color=measurecolors[i])
        if i < length(strs)
            print(io, "\n")
        end
    end
end