
module Associations
    # Use the README as the module docs
    @doc let
        path = joinpath(dirname(@__DIR__), "README.md")
        include_dependency(path)
        read(path, String)
    end Associations

    using Reexport

    using StateSpaceSets
    using DelayEmbeddings: embed, genembed
    export embed, genembed

    import HypothesisTests: pvalue
    export trajectory
    @reexport using StateSpaceSets
    @reexport using ComplexityMeasures
    @reexport using TimeseriesSurrogates

    include("utils/utils.jl")
    include("core.jl")

    include("methods/information/information.jl")
    include("methods/crossmappings/crossmappings.jl")
    include("methods/closeness/closeness.jl")
    include("methods/correlation/correlation.jl")
    include("methods/recurrence/methods.jl")


    # Independence tests must be loaded after everything else has been defined.
    include("independence_tests/independence.jl")

    # Causal graph API must be loaded after independence tests.
    include("causal_graphs/causal_graphs.jl")

    include("deprecations/deprecations.jl")

    # Update messages:
    using Scratch
    display_update = true
    version_number = "4.0.0"
    update_name = "update_v$(version_number)"
    update_message = """
    \nUpdate message: Associations v$(version_number)\n
    - The package has been renamed from CausalityTools.jl to Associations.jl! As part of the renaming, Associations.jl has incremented its major version to v4, which is fully backwards compatible with CausalityTools v3.
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
