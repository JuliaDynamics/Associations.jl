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
end
