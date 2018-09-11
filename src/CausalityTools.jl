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

@reexport using TimeseriesSurrogates
@reexport using StateSpaceReconstruction
@reexport using PerronFrobenius
@reexport using TransferEntropy
using DynamicalSystems

# Example systems
include("systems/Systems.jl")

# Surrogate wrappers for embeddings and Datasets
include("surrogates/surrogates.jl")

end # module
