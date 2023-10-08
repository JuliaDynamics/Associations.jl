export DifferentialDecomposition

# TODO: this should also hold for the dedicated estimators? Or should we have dedicated estimators for those? probably..
"""
    DifferentialDecomposition <: InformationMeasureEstimator
    DifferentialDecomposition(
        definition::MultivariateInformationMeasure,
        est::DifferentialInfoEstimator
    )

The `DifferentialDecomposition` estimates some multivariate information measure by 
decomposing it into a combination of some fundamental/lower-level information measure 
applied to some marginalization of the input data.

For example, [`MIShannon`](@ref) can be decomposed into a sum of entropies, so that 
we can use a [`DifferentialInfoEstimator`](@ref) to 
estimate each of these terms. These estimates can then be combined into a final estimate 
of mutual information. Another example is [`CMIShannon`](@ref), which can both be 
decomposed both into a sum of entropies, and into a sum of mutual informations. This 
means that we can not only  use lower-level entropy estimators
(as for [`MIShannon`](@ref)), but also higher-level [`MutualInformationEstimator`](@ref)s.

!!! note "Bias"
    An estimate obtained through a decomposition into more fundamental information 
    theoretic measures will in general be more biased than when using a dedicated 
    estimator. One reason is that this decomposition may miss out on crucial information
    in the joint space. To remedy this, dedicated estimators typically derives the 
    marginal estimates by first considering the joint space, and then does some clever
    trick to eliminate the bias that is introduced through a naive decomposition.
    For example, [`MutualInformationEstimator`](@ref) are typically more reliable 
    for computing [`MIShannon`](@ref) than using e.g. [`Kraskov`](@ref) as the estimator.

## Controlling the decomposition and estimator

The given estimator `est` controls which decomposition is applied. For example:

- If `est isa DifferentialInfoEstimator`, then 
    the measure is decomposed into its entropy components.
- If `est isa MutualInformationEstimator`, then the measure is decomposed into its
    mutual information components.

!!! info 
    Not all measures have the property that they can be decomposed into more fundamental
    information theoretic quantities. For example, [`MITsallisMartin`](@ref) *can* be 
    decomposed into a combination of marginal entropies, while [`MIRenyiSarbu`](@ref)
    cannot. An error will be thrown if decomposition is not possible.
"""
struct DifferentialDecomposition{M <: MultivariateInformationMeasure, E <: DifferentialInfoEstimator} <: InformationMeasureEstimator{M}
    definition::M # API from complexity measures: definition must be the first field of the info estimator.
    est::E
end

function Base.show(io::IO, est::DifferentialDecomposition)
    types = [
        typeof(est.definition),
        typeof(est.est),
    ]
    strs = [
        "Measure to be decomposed",
        "Estimator for decomposed components",
    ]
    measurecolors = [
        :light_red,
        :light_green,
    ]
    infocolors = [
        :red,
        :green,
    ]
    n = maximum(length.(strs))
    spaces_needed = [n - length(s) for s in strs] 
    spaced_strs = [strs[i] * repeat(" ", spaces_needed[i]) for i in eachindex(strs)]
    ctx = IOContext(io, :color => true)
    printstyled(ctx,  "DifferentialDecomposition estimator\n\n", color=:bold)
    d = decomposition_string(est.definition, est.est)
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
