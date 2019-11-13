module CausalityTools

    using DynamicalSystems
    using StaticArrays

    export Dataset, trajectory, DiscreteDynamicalSystem, ContinuousDynamicalSystem
    export SVector, @SVector, MVector, @MVector
    using Requires

    using Reexport
    @reexport using CausalityToolsBase
    @reexport using StateSpaceReconstruction
    @reexport using PerronFrobenius
    @reexport using TransferEntropy
    @reexport using CrossMappings


    display_update = true
    version = "v0.7.0"
    update_name = "update_$version"

    if display_update
    if !isfile(joinpath(@__DIR__, update_name))
    printstyled(stdout,
    """
    \nUpdate message: CausalityTools $version
    ----------------
    BREAKING CHANGES! 
    ----------------
    
    - The `ν` keyword for the `ConvergentCrossMappingTest` and `CrossMappingTest` has changed to 
        `η`, to conform with the convention for the transfer entropy tests.
    
    There are also some breaking changes to the transfer operator and invariant measure 
    estimators in this release. Other new features:

    - Common interface for causality testing with the `causality` function and its methods.
    - Integration with UncertainData.jl. `causality` accepts uncertain datasets as inputs.
    - Integration with DynamicalSystems.jl. `causality` accepts dynamical systems as inputs.

    Check the online documentation for more information!

    SYNTAX CHANGES
    --------------
      
    Valid inputs
    ------------

    The inputs to the following listed functions should be either a `Dataset` instance, 
    a `CustomReconstruction` instance, or a vector of `Vector`, `SVector` or `MVector`. 

    Invariant measures: 

    - `invariantmeasure(pts, RectangularBinning(ϵ)` where `ϵ` indicates the type of binning.
    - `invariantmeasure(pts, RectangularBinning(ϵ)` where `ϵ` indicates the type of binning.

    Transfer operators:

    - `transferoperator(pts, RectangularBinning(ϵ)` for transfer operators over rectangular 
        binnings.
    - `transferoperator(pts, TriangulationBinning(), ExactIntersection())` for transfer 
        operators over a triangulation of `pts` using exact simplex intersections to compute 
        transition probabilities.
    - `transferoperator(pts, TriangulationBinning(), ApproxIntersection())` for transfer 
    operators over a triangulation of `pts` using approximate simplex intersections to compute 
    transition probabilities.

    """; color = :light_magenta)
    touch(joinpath(@__DIR__, update_name))
    end
    end

    # Example systems
    include("systems/Systems.jl")

    # Wrappers of the different methods.
    include("method_wrappers/highlevel_methods.jl")

    # Various algorithsm that are implemented here and not in subpackages 
    include("algorithms/joint_distance_distribution.jl")
    include("algorithms/smeasure.jl")

    # High-level estimators 
    include("causalitytests/causality_tests.jl")

    ################################################################
    # Integration with DynamicalSystems.jl (for all causality tests)
    ################################################################
    include("integration_dynamicalsystems/IntegrationDynamicalSystems.jl")

    ################################################################
    # Integration with UncertainData.jl (for all causality tests)
    ################################################################
    include("integration_uncertaindata/IntegrationUncertainData.jl")
    
    
    function __init__()
        
        # Surrogate wrappers for embeddings and Datasets
        @require TimeseriesSurrogates="c804724b-8c18-5caa-8579-6025a0767c70" @eval include("surrogates/surrogates.jl")
    end

    # Plot recipes, also for all sub-packages
    #include("plot_recipes/CausalityToolsPlotRecipes.jl")

end # module

