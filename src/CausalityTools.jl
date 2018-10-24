module CausalityTools

using DynamicalSystems
using Distances
using Documenter
using Reexport
import DynamicalSystemsBase:
    DynamicalSystem,
    ContinuousDynamicalSystem,
    DiscreteDynamicalSystem,
    Dataset,
    trajectory

@reexport using TimeseriesSurrogates
@reexport using StateSpaceReconstruction
@reexport using PerronFrobenius
@reexport using TransferEntropy


# Example systems
include("systems/Systems.jl")

# Surrogate wrappers for embeddings and Datasets
include("surrogates/surrogates.jl")

# Wrappers of the different methods.
include("method_wrappers/highlevel_methods.jl")

end # module
