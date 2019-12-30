using Reexport

@reexport module SystemModels

    import DynamicalSystems: ContinuousDynamicalSystem, AbstractDataset, Dataset, trajectory
    import StaticArrays: SVector 
    import StatsBase: std 
    import Distributions: Normal 

    # Abstract model definitions and methods
    include("AbstractSystemModel.jl")
    include("ContinuousSystemModel.jl")
    include("DiscreteSystemModel.jl")

    # Continuous system definitions
    include("continuous/RosslerLorenzUnidir.jl")

    export 
    AbstractSystemModel, 
    ContinuousSystemModel,
    DiscreteSystemModel,
    randomised,
    get_dt,
    get_ui,
    get_nvars,
    trajectory,
    RosslerLorenzUnidir

end