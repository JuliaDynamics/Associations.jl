module CausalityTools

    using DynamicalSystems
    using StaticArrays

    export Dataset, trajectory, DiscreteDynamicalSystem, ContinuousDynamicalSystem
    export SVector, @SVector, MVector, @MVector

    using Reexport
    @reexport using TimeseriesSurrogates
    @reexport using CausalityToolsBase
    @reexport using StateSpaceReconstruction
    @reexport using PerronFrobenius
    @reexport using TransferEntropy
    @reexport using CrossMappings


    display_update = true
    version = "v0.3.0"
    update_name = "update_$version"

    if display_update
    if !isfile(joinpath(@__DIR__, update_name))
    printstyled(stdout,
    """
    \nUpdate message: CausalityTools $version
    ----------------
    BREAKING CHANGES! 
    ----------------

    VALID INPUTS 
    The inputs to the following listed functions should be either a `Dataset` instance, 
    a `CustomReconstruction` instance, or a vector of `Vector`, `SVector` or `MVector`. 

    SYNTAX CHANGES
    --------------
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

    # Surrogate wrappers for embeddings and Datasets
    include("surrogates/surrogates.jl")

    # Wrappers of the different methods.
    include("method_wrappers/highlevel_methods.jl")

    # Various algorithsm that are implemented here and not in subpackages 
    include("algorithms/joint_distance_distribution.jl")
    
    # High-level estimators 
    include("causalitytests/causalitytests.jl")

    ################################################################
    # Integration with DynamicalSystems.jl (for all causality tests)
    ################################################################
    include("integration_dynamicalsystems/IntegrationDynamicalSystems.jl")

    ################################################################
    # Integration with UncertainData.jl (for all causality tests)
    ################################################################
    include("integration_uncertaindata/IntegrationUncertainData.jl")
    
    # Plot recipes, also for all sub-packages
    #include("plot_recipes/CausalityToolsPlotRecipes.jl")

end # module

