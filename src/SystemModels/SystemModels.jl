using Reexport

@reexport module SystemModels

    import DynamicalSystems: 
        ContinuousDynamicalSystem, 
        AbstractDataset, 
        Dataset, 
        trajectory,
        columns
    using UncertainData
    import StaticArrays: SVector 
    import StatsBase: std 
    import SimpleDiffEq
    import SimpleDiffEq: get_dt
    using Distributions
    import Base.rand

    # Abstract model definitions and methods
    include("AbstractSystemModel.jl")
    include("ContinuousSystemModel.jl")
    include("DiscreteSystemModel.jl")

    # Methods to add observational noise after orbits have been obtained 
    include("add_noise.jl")

    # Continuous system definitions
    include("continuous/RosslerLorenzUnidir.jl")

    export 
    AbstractSystemModel, 
    ContinuousSystemModel,
    DiscreteSystemModel,
    rand,
    get_dt,
    get_ui,
    get_nvars,
    add_observational_noise!,
    trajectory,
    RosslerLorenzUnidir

end