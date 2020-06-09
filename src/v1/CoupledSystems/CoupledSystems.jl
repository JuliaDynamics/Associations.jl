using Reexport 

@reexport module CoupledSystems
    using StaticArrays
    using DynamicalSystemsBase

    include("continuous/lorenzdiffusive.jl")
    include("discrete/latticemapunidir.jl")
end