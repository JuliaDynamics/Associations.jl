using CausalityTools
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

using DynamicalSystems
import DynamicalSystemsBase:
    DynamicalSystem,
    DiscreteDynamicalSystem,
    ContinuousDynamicalSystem,
    Dataset

@testset "Discrete systems" begin
	include("discrete_systems.jl")
end

@testset "Continuous systems" begin
	include("continuous_systems.jl")
end

@testset "Transfer entropy wrappers" begin
	include("wrappers_te.jl")
end
