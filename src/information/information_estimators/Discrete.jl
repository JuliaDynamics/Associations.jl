export DiscreteDecomposition
export codified_marginals

"""
    DiscreteDecomposition <: InformationMeasureEstimator
    DiscreteDecomposition(definition::MultivariateInformationMeasure,
        est::DiscreteInfoEstimator, 
        discretization, 
        pest::ProbabilitiesEstimator = RelativeAmount())

Like [`DifferentialDecomposition`](@ref), but uses discrete estimation. 

The given `discretization` scheme (typically an [`OutcomeSpace`](@ref)) controls how the
joint/marginals are discretized, and the probabilities estimator `pest` controls how
probabilities are estimated from counts.

!!! note "Bias"
    Like for [`DifferentialDecomposition`](@ref), using a dedicated estimator 
    for the measure in question will be more reliable than using a decomposition
    estimate. Here's how different `discretization`s are applied:

    - [`ValueBinning`](@ref). Bin visitation frequencies are counted in the joint space `XY`,
        then marginal visitations are obtained from the joint bin visits.
        This behaviour is the same for both [`FixedRectangularBinning`](@ref) and
        [`RectangularBinning`](@ref) (which adapts the grid to the data).
        When using [`FixedRectangularBinning`](@ref), the range along the first dimension
        is used as a template for all other dimensions. This is a bit slower than naively 
        binning each marginal, but lessens bias.
    - [`OrdinalPatterns`](@ref). Each timeseries is separately [`codify`](@ref)-ed according
        to its ordinal pattern (so bias is presumably high).
    - [`Dispersion`](@ref). Each timeseries is separately [`codify`](@ref)-ed according to its
        dispersion pattern  (so bias is presumably high).

## Examples

```julia
using CausalityTools
using Random; rng = MersenneTwister(1234)

x = StateSpaceSet(rand(rng, 1000000, 2))
y = StateSpaceSet(rand(rng, 1000000, 2))
# Compute Shannon mutual information by discretizing each marginal column-wise
# (per variable) using length-`3` ordinal patterns.
est = DiscreteDecomposition(MIShannon(), PlugIn(Shannon()), OrdinalPatterns(m=2))
information(est, x, y) # should be close to 0
```
"""
struct DiscreteDecomposition{M <: MultivariateInformationMeasure, E, D, P} <: InformationMeasureEstimator{M}
    definition::M # API from complexity measures: definition must be the first field of the info estimator.
    mest::E # The estimator + measure which `est.definition` is decomposed into.
    discretization::D
    pest::P # How probabilities are computed

    function DiscreteDecomposition(definition::M, mest::E, discretization::D, pest::P = RelativeAmount()) where {M <: MultivariateInformationMeasure, E, D, P}
        return new{M, E, D, P}(definition, mest, discretization, pest)
    end
end

function decomposition_string(
        definition::MultivariateInformationMeasure, 
        est::InformationMeasureEstimator,
    )
    return "The formula for decomposition of $definition is not available. " * 
        "Please open an issue.";
end

function Base.show(io::IO, est::DiscreteDecomposition)

    types = [
        typeof(est.definition),
        typeof(est.mest),
        typeof(est.discretization),
        typeof(est.pest),
    ]
    strs = [
        "Measure to be decomposed",
        "Estimator for decomposed components",
        "Discretization",
        "Probabilities estimator",
    ]
    measurecolors = [
        :light_red,
        :light_green,
        :light_yellow,
        :light_blue
    ]
    infocolors = [
        :red,
        :green,
        :yellow,
        :blue
    ]
    n = maximum(length.(strs))
    spaces_needed = [n - length(s) for s in strs] 
    spaced_strs = [strs[i] * repeat(" ", spaces_needed[i]) for i in eachindex(strs)]
    ctx = IOContext(io, :color => true)
    printstyled(ctx,  "DiscreteDecomposition estimator\n\n", color=:bold)
    d = decomposition_string(est.definition, est.mest)
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


"""
    codified_marginals(o::OutcomeSpace, x::VectorOrStateSpaceSet...)

Encode/discretize each input vector `xᵢ ∈ x` according to a procedure determined by `o`.
Any `xᵢ ∈ X` that are multidimensional ([`StateSpaceSet`](@ref)s) will be encoded column-wise,
i.e. each column of `xᵢ` is treated as a timeseries and is encoded separately.

This is useful for computing any discrete information theoretic quantity.

## Supported estimators

- [`ValueBinning`](@ref). Bin visitation frequencies are counted in the joint space `XY`,
    then marginal visitations are obtained from the joint bin visits.
    This behaviour is the same for both [`FixedRectangularBinning`](@ref) and
    [`RectangularBinning`](@ref) (which adapts the grid to the data).
    When using [`FixedRectangularBinning`](@ref), the range along the first dimension
    is used as a template for all other dimensions.
- [`OrdinalPatterns`](@ref). Each timeseries is separately [`encode`](@ref)d according
    to its ordinal pattern.
- [`Dispersion`](@ref). Each timeseries is separately [`encode`](@ref)d according to its
    dispersion pattern.

Many more implementations are possible. Each new implementation gives one new
way of estimating the [`ContingencyMatrix`](@ref)
"""
function codified_marginals end

function codified_marginals(est, x::VectorOrStateSpaceSet...)
    return codify_marginal.(Ref(est), x)
end

function codify_marginal(est, x::AbstractStateSpaceSet)
    return StateSpaceSet(codify_marginal.(Ref(est), columns(x))...)
end

function codify_marginal(o::UniqueElements, x::AbstractVector)
    return x
end

function codify_marginal(o::OrdinalPatterns{m}, x::AbstractVector) where {m}
    return codify(o, x)
end

function codify_marginal(o::Dispersion, x::AbstractVector)
    return codify(o, x)
end

function codify_marginal(
        o::ValueBinning{<:FixedRectangularBinning{D}},
        x::AbstractVector) where D
    range = first(o.binning.ranges)
    ϵmin = minimum(range)
    ϵmax = maximum(range)
    N = length(range)
    encoder = RectangularBinEncoding(FixedRectangularBinning(ϵmin, ϵmax, N, 1))
    return encode.(Ref(encoder), x)
end

# Special treatment for RectangularBinning. We create the joint embedding, then
# extract marginals from that. This could probably be faster,
# but it *works*. I'd rather things be a bit slower than having marginals
# that are not derived from the same joint distribution, which would hugely increase
# bias, because we're not guaranteed cancellation between entropy terms
# in higher-level methods.
function codified_marginals(o::ValueBinning{<:RectangularBinning}, x::VectorOrStateSpaceSet...)
    # TODO: The following line can be faster by explicitly writing out loops that create the 
    # joint embedding vectors.
    X = StateSpaceSet(StateSpaceSet.(x)...)
    encoder = RectangularBinEncoding(o.binning, X)

    bins = [vec(encode_as_tuple(encoder, pt))' for pt in X]
    joint_bins = reduce(vcat, bins)
    idxs = size.(x, 2) #each input can have different dimensions
    s = 1
    encodings = Vector{Vector}(undef, 0)
    for (i, cidx) in enumerate(idxs)
        variable_subset = s:(s + cidx - 1)
        s += cidx
        y = @views joint_bins[:, variable_subset]
        for j in size(y, 2)
            push!(encodings, y[:, j])
        end
    end

    return encodings
end

# A version of `cartesian_bin_index` that directly returns the joint bin encoding
# instead of converting it to a cartesian index.
function encode_as_tuple(e::RectangularBinEncoding, point::SVector{D, T}) where {D, T}
    ranges = e.ranges
    if e.precise
        # Don't know how to make this faster unfurtunately...
        bin = map(searchsortedlast, ranges, point)
    else
        bin = floor.(Int, (point .- e.mini) ./ e.widths) .+ 1
    end
    return bin
end