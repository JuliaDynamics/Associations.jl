module CausalityTools

using DynamicalSystems
using Documenter
using Reexport
import DynamicalSystemsBase:
    DynamicalSystem,
    ContinuousDynamicalSystem,
    DiscreteDynamicalSystem,
    Dataset,
    trajectory

@reexport using StateSpaceReconstruction
@reexport using PerronFrobenius
@reexport using TransferEntropy
@reexport using DynamicalSystems.trajectory
@reexport using DynamicalSystems.Dataset

# Example systems
include("systems/Systems.jl")

end # module
