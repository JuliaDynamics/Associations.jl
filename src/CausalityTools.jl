module CausalityTools
    using Reexport
    import DelayEmbeddings: Dataset
    import DynamicalSystemsBase: trajectory
    import DynamicalSystemsBase: DiscreteDynamicalSystem, ContinuousDynamicalSystem
    export Dataset
    export trajectory
    export DiscreteDynamicalSystem, ContinuousDynamicalSystem
    @reexport using Entropies
    @reexport using TimeseriesSurrogates

    include("core.jl")
    include("methods/infomeasures/infomeasures.jl")
    include("methods/ccm/ccm.jl")
    include("methods/joint_distance_distribution/joint_distance_distribution.jl")
    include("methods/pairwise_asymmetric_inference/pairwise_asymmetric_inference.jl")
    include("methods/smeasure/smeasure.jl")
    include("example_systems/ExampleSystems.jl")

    include("utils/kde.jl")

    # Independence tests must be loaded after everything else has been defined.
    include("independence.jl")

    #using Requires
    #function __init__()
        #@require UncertainData="dcd9ba68-c27b-5cea-ae21-829cd07325bf" begin
        #   include("integrations/uncertaindata.jl")
        #end

        #@require Simplices="d5428e67-3037-59ba-9ab1-57a04f0a3b6a" begin
        # #   import PerronFrobenius: SimplexExact, SimplexPoint
        #    export SimplexExact, SimplexPoint
        #end
    #end

    # Update messages:
    using Scratch
    display_update = true
    version_number = "2.0.0"
    update_name = "update_v$(version_number)"
    update_message = """
    \nUpdate message: CausalityTools v$(version_number)\n
    - An overall overhaul of the documentation and API of CausalityTools.jl has been
        performed.
    - The syntax for all information-based methods have changed. Estimators, which
        also contains analysis parameters, are now always the first argument.
    - All information-based methods in the DynamicalSystems.jl organization that are
        more complex than those in `Entropies.jl` have been moved to CausalityTools.jl.
    - A lot of new methods and estimators have been added. See the online documentation for
       an overview.
    - Functionality from TransferEntropy.jl has been moved to CausalityTools.jl.
    - The `Kraskov1` and `Kraskov2` MI estimators have been renamed to
        `KraskovStögbauerGrassberger1` (`KSG1` for short) and
        `KraskovStögbauerGrassberger2` (`KSG2` for short),
        and can now also compute MI between more than two datasets at once.
    """

    if display_update
        # Get scratch space for this package
        versions_dir = @get_scratch!("versions")
        if !isfile(joinpath(versions_dir, update_name))
            printstyled(
                stdout,
                update_message;
                color = :light_magenta,
            )
            touch(joinpath(versions_dir, update_name))
        end
    end
end
