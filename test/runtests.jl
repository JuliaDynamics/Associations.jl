if lowercase(get(ENV, "CI", "false")) == "true"
    include("install_dependencies.jl")
end

using Test
using DynamicalSystems
import DynamicalSystemsBase:
    DynamicalSystem,
    DiscreteDynamicalSystem,
    ContinuousDynamicalSystem,
    Dataset

using StateSpaceReconstruction
using PerronFrobenius
using CausalityTools

@testset "Plot recipes" begin
	include("plot_recipes.jl")
end

@testset "Discrete systems" begin
	include("discrete_systems.jl")
end

@testset "Continuous systems" begin
	include("continuous_systems.jl")
end

@testset "Transfer entropy wrappers" begin
	include("wrappers_te.jl")
end
