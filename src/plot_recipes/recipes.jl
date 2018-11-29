# Helper functions used in the various recipes. These must be imported first.
include("helperfunctions_rectangulargrid.jl")
include("helperfunctions_triangulationgrid.jl")

# The actual recipes.
using StateSpaceReconstruction
include("recipe_embedding.jl")

using PerronFrobenius
include("recipe_rectangularinvariantmeasure.jl")
include("recipe_invariantdistribution.jl")
include("recipe_transferoperator.jl")
include("recipe_triangulationplot.jl")
