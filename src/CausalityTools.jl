
module CausalityTools
    # Use the README as the module docs
    @doc let
        path = joinpath(dirname(@__DIR__), "README.md")
        include_dependency(path)
        read(path, String)
    end CausalityTools

    using Reexport

    using StateSpaceSets
    using DelayEmbeddings: embed, genembed
    export embed, genembed

    import DynamicalSystemsBase: trajectory
    import DynamicalSystemsBase: DiscreteDynamicalSystem, ContinuousDynamicalSystem
    import HypothesisTests: pvalue
    export trajectory
    export DiscreteDynamicalSystem, ContinuousDynamicalSystem
    @reexport using StateSpaceSets
    @reexport using ComplexityMeasures
    @reexport using TimeseriesSurrogates

    include("core.jl")
    include("methods/infomeasures/infomeasures.jl")
    include("methods/crossmappings/crossmappings.jl")
    include("methods/closeness/closeness.jl")
    include("methods/correlation/correlation.jl")
    include("methods/recurrence/methods.jl")

    include("utils/utils.jl")

    # Independence tests must be loaded after everything else has been defined.
    include("independence_tests/independence.jl")

    # Causal graph API must be loaded after independence tests.
    include("causal_graphs/causal_graphs.jl")

    include("example_systems/example_systems.jl")

    include("deprecations/deprecations.jl")

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
    version_number = "2.6.0"
    update_name = "update_v$(version_number)"
    update_message = """
    \nUpdate message: CausalityTools v$(version_number)\n
    - New causal graph inference algorithm: `PC`.
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
