#if lowercase(get(ENV, "CI", "false")) == "true"
#    include("install_dependencies.jl")
#end

using Test
using CausalityTools

@testset "SMeasure" begin 
    x, y = rand(100), rand(100)
    @test s_measure(x, y) isa Float64
end

@testset "PredictiveAsymmetry" begin 
    x, y, z = rand(100), rand(100), rand(100)

    # Define prediction lags and estimation method
    ηs = 1:5
    method = VisitationFrequency(RectangularBinning(5))

    # 𝔸(x → y) and  𝔸(x → y | z)
    𝔸reg  = predictive_asymmetry(x, y, method, ηs, normalize = false)
    𝔸cond = predictive_asymmetry(x, y, z, method, ηs, normalize = false)
    𝒜reg = predictive_asymmetry(x, y, method, ηs,  f = 1.0) # normalize == true by default
    𝒜cond = predictive_asymmetry(x, y, z, method, ηs, f = 1.5) # normalize == true by default

    @test 𝔸reg isa Vector{<:Real}
    @test 𝔸cond isa Vector{<:Real}
    @test 𝒜reg isa Vector{<:Real}
    @test 𝒜cond isa Vector{<:Real}
end

@testset "Discrete example systems" begin 
    include("systems/discrete/test_discrete_systems.jl")
end

@testset "Continuous example systems" begin 
    include("systems/continuous/test_continuous_systems.jl")
end

#using Distributions
#using UncertainData#

# @testset "Plot recipes" begin
# 	include("plot_recipes.jl")
# end

# @testset "SystemModels" begin
#     include("SystemModels/test_SystemModels.jl")
# end

# @testset "CausalAnalysis" begin
#     include("causality_tests/CausalAnalysis/test_CausalAnalysis.jl")
#     include("causality_tests/CausalAnalysis/test_CausalAnalysis_meta.jl")
# end



# @testset "High level wrappers" begin
#     include("test_wrappers_te.jl")
# end

# @testset "Causality tests on scalar time series" begin
#     include("causality_tests/test_scalarseries_CrossMappingTest.jl")
#     include("causality_tests/test_scalarseries_ConvergentCrossMappingTest.jl")
#     include("causality_tests/test_scalarseries_JointDistanceDistributionTest.jl")
#     include("causality_tests/test_scalarseries_JointDistanceDistributionTTest.jl")
#     include("causality_tests/test_scalarseries_TransferOperatorGridTest.jl")
#     include("causality_tests/test_scalarseries_VisitationFrequencyTest.jl")
#     include("causality_tests/test_scalarseries_NearestNeighbourMITest.jl")

#     include("causality_tests/test_scalarseries_ApproximateSimplexIntersectionTest.jl")
#     include("causality_tests/test_scalarseries_ExactSimplexIntersectionTest.jl")
#     include("causality_tests/test_scalarseries_PredictiveAsymmetryTest.jl")
#     include("causality_tests/test_scalarseries_NormalisedPredictiveAsymmetryTest.jl")
#     include("causality_tests/test_scalarseries_SMeasure.jl")
# end


# @testset "DynamicalSystems.jl integration" begin
#     include("causality_tests/integration_dynamicalsystems/test_integration_discrete_system.jl")
#     include("causality_tests/integration_dynamicalsystems/test_integration_continuous_system.jl")
# end

# @testset "UncertainData.jl integration" begin
#     include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_JointDistanceDistributionTest.jl")
#     include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_JointDistanceDistributionTTest.jl")
#     include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_CrossMappingTest.jl")
#     include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_ConvergentCrossMappingTest.jl")
#     include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_TransferOperatorGridTest.jl")
#     include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_VisitationFrequencyTest.jl")
#     include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_ApproximateSimplexIntersectionTest.jl")
#     include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_ExactSimplexIntersectionTest.jl")
#     include("causality_tests/integration_uncertaindata/test_uncertaindata_integration_PredictiveAsymmetryTest.jl")

#     include("causality_tests/integration_uncertaindata/test_UncertainIndexValueDataset.jl")
    
#     # High-level tests
#     include("causality_tests/highlevel_tests/test_ConstrainedTest.jl")
#     include("causality_tests/highlevel_tests/test_InterpolateBinTest.jl")
#     include("causality_tests/highlevel_tests/test_RandomSequencesTest.jl")

#     include("causality_tests/integration_uncertaindata/test_uncertain_indexvalue_dataset_with_schemes.jl")

# end
