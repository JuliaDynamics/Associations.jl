module CausalityTools

using DynamicalSystems
using Distances
using Reexport
using AbstractFFTs
using Statistics
using RecipesBase
using StaticArrays
using LinearAlgebra
using Measures


@reexport using TimeseriesSurrogates
@reexport using StateSpaceReconstruction
@reexport using PerronFrobenius
@reexport using TransferEntropy
@reexport using CrossMappings

using StateSpaceReconstruction: GroupSlices

import DynamicalSystemsBase:
    DynamicalSystem,
    ContinuousDynamicalSystem,
    DiscreteDynamicalSystem,
    Dataset,
    trajectory

# Example systems
include("systems/Systems.jl")

# Surrogate wrappers for embeddings and Datasets
include("surrogates/surrogates.jl")

# Wrappers of the different methods.
include("method_wrappers/highlevel_methods.jl")

# Plot recipes, also for all sub-packages
include("plot_recipes/recipes.jl")

export
Dataset,
DynamicalSystem,
DiscreteDynamicalSystem,
ContinuousDynamicalSystem,
trajectory

end # module
