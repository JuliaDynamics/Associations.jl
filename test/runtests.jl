if lowercase(get(ENV, "CI", "false")) == "true"
    include("install_dependencies.jl")
end

using Test
using CausalityTools
using Distributions
using UncertainData

#@testset "Plot recipes" begin
#	include("plot_recipes.jl")
#end

@testset "Discrete example systems" begin 
    include("systems/discrete/test_discrete_systems.jl")
end

@testset "Continuous example systems" begin 
    include("systems/continuous/test_continuous_systems.jl")
end

@testset "High level wrappers" begin
    include("test_wrappers_te.jl")
end

@testset "Causality tests on scalar time series" begin
    include("causality_tests/test_scalarseries_CrossMappingTest.jl")
    include("causality_tests/test_scalarseries_ConvergentCrossMappingTest.jl")
    include("causality_tests/test_scalarseries_JointDistanceDistributionTest.jl")
    include("causality_tests/test_scalarseries_JointDistanceDistributionTTest.jl")
    include("causality_tests/test_scalarseries_TransferOperatorGridTest.jl")
    include("causality_tests/test_scalarseries_VisitationFrequencyTest.jl")
    include("causality_tests/test_scalarseries_ApproximateSimplexIntersectionTest.jl")
    include("causality_tests/test_scalarseries_PredictiveAsymmetryTest.jl")
end

@testset "DynamicalSystems.jl integration" begin
    include("causality_tests/integration_dynamicalsystems/test_integration_discrete_system.jl")
    include("causality_tests/integration_dynamicalsystems/test_integration_continuous_system.jl")
end

@testset "UncertainData.jl integration" begin
    include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_JointDistanceDistributionTest.jl")
    include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_JointDistanceDistributionTTest.jl")
    include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_CrossMappingTest.jl")
    include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_ConvergentCrossMappingTest.jl")
    include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_TransferOperatorGridTest.jl")
    include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_VisitationFrequencyTest.jl")
    include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_ApproximateSimplexIntersectionTest.jl")
    include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_PredictiveAsymmetryTest.jl")
end
