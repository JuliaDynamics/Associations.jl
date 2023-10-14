# Exists just to make dispatch easier for pretty printing.
abstract type DecompositionEstimator{M} <: MultivariateInformationMeasureEstimator{M} end

include("EntropyDecomposition.jl")
include("MIDecomposition.jl")
include("CMIDecomposition.jl")

# ----------------------------------------------------------------------------------------
# Pretty printing
# ----------------------------------------------------------------------------------------
# A common method for displaying. For pretty printing to to, each `DecompositionEstimator`
# must implement relevant methods for `summary_types`, `summary_strings`,
# `measure_colors` and `info_colors`.
# 
# If custom printing isn't defined for a particular measure-definition combination,
# then default type printing is used.
function Base.show(io::IO, est::DecompositionEstimator)
    types = summary_types(est)
    strs = summary_strings(est)
    measurecolors = measure_colors(est)
    infocolors = info_colors(est)
    n = maximum(length.(strs))

    spaces_needed = [n - length(s) for s in strs] 
    spaced_strs = [strs[i] * repeat(" ", spaces_needed[i]) for i in eachindex(strs)]
    ctx = IOContext(io, :color => true)
    printstyled(ctx,  "$(typeof(est).name.name) estimator\n\n", color=:bold)
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

function summary_strings(est::DecompositionEstimator)
    return [
        "Measure to be decomposed",
        "Estimator for components",
    ]
end

function summary_types(est::DecompositionEstimator)
    return [
        typeof(est.definition),
        typeof(est.est),
    ]
end

function measure_colors(est::DecompositionEstimator)
    return [
        :light_red,
        :light_green,
    ]
end

function info_colors(est::DecompositionEstimator)
    return [
        :red,
        :green,
    ]
end
