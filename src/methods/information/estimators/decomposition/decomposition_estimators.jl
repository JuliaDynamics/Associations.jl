# Exists just to make dispatch easier for pretty printing.
abstract type DecompositionEstimator{M} <: MultivariateInformationMeasureEstimator{M} end

include("EntropyDecomposition.jl")
include("MIDecomposition.jl")
include("CMIDecomposition.jl")

function decomposition_string(::SOURCE_DEF, ::DecompositionEstimator{<:TARGET_DEF}) where {SOURCE_DEF, TARGET_DEF}
    #"Decomposition formula for $M not defined for $E"
    "Not specified. Are you sure $SOURCE_DEF is the measure you want to decompose " * 
    "into $TARGET_DEF?"
end
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

    max_str_width = maximum(length.(strs))
    max_type_width = maximum(length.(string.(types)))
    formula = decomposition_string(est.definition, est)

    max_content_width = maximum([
        length("$(typeof(est).name.name) estimator"),
        max_str_width + max_type_width + 2, 
        length(formula) + 10  # +10 for " Formula: "
    ])

    ctx = IOContext(io, :color => true)

    function print_horizontal_line(char::Char)
        printstyled(ctx, "+$(repeat(char, max_content_width + 2))+\n", color = :bold)
    end

    print_horizontal_line('-')
    printstyled(ctx, "| $(center_str("$(typeof(est).name.name) estimator", max_content_width)) |\n", color = :bold)
    print_horizontal_line('-')
   
    for i in eachindex(strs)
        printstyled(ctx, "| ", color = :bold)
        printstyled(ctx, "$(rpad(strs[i], max_str_width)): ", color = :grey)
        type_str = string(types[i])
        printstyled(ctx, type_str, color = measurecolors[i])
        padding = max_content_width - max_str_width - length(type_str) - 2
        printstyled(ctx, repeat(" ", padding), color = :bold)
        printstyled(ctx, " |\n", color = :bold)
    end

    formula_line = "Formula: $(formula)"
    printstyled(ctx, "| ", color = :bold)
    printstyled(ctx, formula_line, color = :light_grey)
    padding = max_content_width - length(formula_line)
    printstyled(ctx, repeat(" ", padding), color = :bold)
    printstyled(ctx, " |\n", color = :bold)
    print_horizontal_line('-')
end

function center_str(s::String, width::Int)
    if length(s) >= width
        return s
    else
        lpad = (width - length(s)) ÷ 2
        rpad = width - length(s) - lpad
        return "$(repeat(" ", lpad))$(s)$(repeat(" ", rpad))"
    end
end

function summary_strings(est::DecompositionEstimator)
    return [
        "Measure",
        "Component estimator",
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
