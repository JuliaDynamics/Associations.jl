module CausalityTools
    using Reexport 
    import DynamicalSystemsBase: trajectory, DiscreteDynamicalSystem, ContinuousDynamicalSystem
    import DelayEmbeddings: Dataset
    export trajectory, DiscreteDynamicalSystem, ContinuousDynamicalSystem, Dataset
    
    import TransferEntropy
    import TransferEntropy: transferentropy, mutualinfo, Hilbert
    export transferentropy, mutualinfo, Hilbert
    @reexport using Entropies
    @reexport using TransferEntropy
    @reexport using TimeseriesSurrogates

    include("JointDistanceDistribution/JointDistanceDistribution.jl")
    include("CrossMappings/CrossMappings.jl")
    include("SMeasure/smeasure.jl")
    include("PredictiveAsymmetry/PredictiveAsymmetry.jl")
    include("Leanings/Leanings.jl")
    include("example_systems/ExampleSystems.jl")

    using Requires 
    function __init__()
        #@require UncertainData="dcd9ba68-c27b-5cea-ae21-829cd07325bf" begin
        #   include("integrations/uncertaindata.jl")
        #end

        #@require Simplices="d5428e67-3037-59ba-9ab1-57a04f0a3b6a" begin
        # #   import PerronFrobenius: SimplexExact, SimplexPoint
        #    export SimplexExact, SimplexPoint
        #end
    end
end 

