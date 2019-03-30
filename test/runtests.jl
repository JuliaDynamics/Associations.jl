if lowercase(get(ENV, "CI", "false")) == "true"
    include("install_dependencies.jl")
end

using Test
using CausalityTools

#@testset "Plot recipes" begin
#	include("plot_recipes.jl")
#end

include("test_discrete_systems.jl")
include("test_continuous_systems.jl")

@testset "High level wrappers" begin
	include("test_wrappers_te.jl")
end
