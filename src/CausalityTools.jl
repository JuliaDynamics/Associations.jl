module CausalityTools
    using Reexport 
    import DynamicalSystemsBase: trajectory, DiscreteDynamicalSystem, ContinuousDynamicalSystem
    export trajectory, DiscreteDynamicalSystem, ContinuousDynamicalSystem
    
    import TransferEntropy
    import TransferEntropy: transferentropy, mutualinfo, Hilbert
    export transferentropy, mutualinfo, Hilbert
    @reexport using Entropies
    @reexport using TransferEntropy
    @reexport using TimeseriesSurrogates
    @reexport using CrossMappings

    include("JointDistanceDistribution/JointDistanceDistribution.jl")
    include("SMeasure/smeasure.jl")
    include("PredictiveAsymmetry/PredictiveAsymmetry.jl")
    include("example_systems/ExampleSystems.jl")

    using Requires 
    function __init__()
        @require UncertainData="dcd9ba68-c27b-5cea-ae21-829cd07325bf" begin
            include("integrations/uncertaindata.jl")
        end

        #@require Simplices="d5428e67-3037-59ba-9ab1-57a04f0a3b6a" begin
        # #   import PerronFrobenius: SimplexExact, SimplexPoint
        #    export SimplexExact, SimplexPoint
        #end
    end
    #using DynamicalSystems
    #using StaticArrays

    #export Dataset, trajectory, DiscreteDynamicalSystem, ContinuousDynamicalSystem
    #export SVector, @SVector, MVector, @MVector
    #using Requires

    #using Reexport
    #@reexport using CausalityToolsBase
    #@reexport using StateSpaceReconstruction
    #@reexport using PerronFrobenius
    #@reexport using TransferEntropy
    #@reexport using CrossMappings

    # Modules defined in CausalityTools
    # include("JointDistanceDistribution/JointDistanceDistribution.jl")

    # display_update = true
    # version = "v0.9.0"
    # update_name = "update_$version"

    # if display_update
    # if !isfile(joinpath(@__DIR__, update_name))
    # printstyled(stdout,
    # """
    # \nUpdate message: CausalityTools $version
    # ----------------
    # New functionality  
    # ----------------
    # #### TransferEntropy.jl

    # - Added [`NearestNeighbourMI`](@ref) transfer entropy estimator.
    # - Added [`transferentropy(::Any, ::TEVars, ::NearestNeighbourMI)`](@ref) method.
    # - Transfer entropy estimators now contain a field `b` which gives the base of the logarithm
    #     used during transfer entropy computations, and hence dictates the unit of the transfer 
    #     entropy. By default, `b = 2`, which gives the transfer entropy in bits. The keyword `b` 
    #     is thus obsolete in all transfer entropy methods that used it before.

    # #### CausalityTools.jl

    # - Added [`NearestNeighbourMITest`](@ref) transfer entropy test.
    # - Added [`NormalisedPredictiveAsymmetryTest`](@ref), for which the predictive asymmetry 
    #     is normalised to some fraction of the mean of the raw values of the predictive statistic.
    # - All subtypes of `CausalityTest` are now mutable. This allows adjusting test 
    #     parameters during sensitivity analyses without creating a new test every time.

    # #### TimeseriesSurrogates.jl

    # - Plotting functionality is behind a Requires-block. To use the plotting functionality, you 
    #     now need to do `using Plots` beforehand.

    # ### Documentation

    # - Fixed small error in documentation for `DiscreteSystemSetup`.
    # - Created high-level pages for transfer entropy. Re-did documentation.
    # - Created high-level pages for cross mappings.
    # """; color = :light_magenta)
    # touch(joinpath(@__DIR__, update_name))
    # end
    # end

    # # Example systems
    # include("systems/Systems.jl")

    # # SystemModels module
    # include("SystemModels/SystemModels.jl")

    # # Wrappers of the different methods.
    # include("method_wrappers/highlevel_methods.jl")
    # include("method_wrappers/transferentropy_from_timeseries.jl")

    # # Various algorithsm that are implemented here and not in subpackages 
    # include("algorithms/smeasure.jl")

    # # High-level estimators 
    # include("CausalityTests/CausalityTests.jl")

    # ################################################################
    # # Integration with DynamicalSystems.jl (for all causality tests)
    # ################################################################
    # include("integration_dynamicalsystems/IntegrationDynamicalSystems.jl")

    # ################################################################
    # # Integration with UncertainData.jl (for all causality tests)
    # ################################################################
    # include("integration_uncertaindata/IntegrationUncertainData.jl")
    
    
    # function __init__()
        
    #     # Surrogate wrappers for embeddings and Datasets
    #     @require TimeseriesSurrogates="c804724b-8c18-5caa-8579-6025a0767c70" @eval include("surrogates/surrogates.jl")
    # end

    # # Plot recipes, also for all sub-packages
    # #include("plot_recipes/CausalityToolsPlotRecipes.jl")

end # module

