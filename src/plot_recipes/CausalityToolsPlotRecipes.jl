using Reexport

@reexport module CausalityToolsPlotRecipes
    using RecipesBase

    # Helper functions used in the various recipes. These must be imported first.
    include("helperfunctions_rectangulargrid.jl")
    include("helperfunctions_triangulationgrid.jl")

    # The actual recipes.
    include("recipe_embedding.jl")
    include("recipe_rectangularinvariantmeasure.jl")
    include("recipe_invariantdistribution.jl")
    include("recipe_transferoperator.jl")
    include("recipe_Simplices.jl")
    include("recipe_triangulationplot.jl")

end # module
