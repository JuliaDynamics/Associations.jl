using Reexport

@reexport module EmpiricalDynamicalModelling

    include("simplex_projection.jl")
    include("delay_simplexprojection.jl")
    include("smap.jl")
end