using CausalityTools
using RecipesBase


@testset "Recipes for subpackages" begin

    # Test plot recipes by calling RecipesBase.apply_recipe with empty dict.
    # It should return a vector of RecipesBase.RecipeData
    d = Dict{Symbol,Any}()

    # Create an example embedding
    E = StateSpaceReconstruction.embed([diff(rand(30)) for i = 1:3], [1, 2, 3], [1, 0, -1])
    E_invariant = invariantize(E)
    invm = RectangularInvariantMeasure(E, [4, 5, 6])


    @testset "StateSpaceReconstruction.AbstractEmbedding" begin
        @test typeof(RecipesBase.apply_recipe(d, E)) == Array{RecipesBase.RecipeData,1}
        @test typeof(RecipesBase.apply_recipe(d, E_invariant)) == Array{RecipesBase.RecipeData,1}
    end

    @testset "PerronFrobenius." begin
        typeof(RecipesBase.apply_recipe(d, invm.transfermatrix)) == Array{RecipesBase.RecipeData,1}
    end

    @testset "Invariant distribution" begin
        typeof(RecipesBase.apply_recipe(d, invm.measure)) == Array{RecipesBase.RecipeData,1}
    end

    @testset "RectangularInvariantMeasure" begin
        boxfillfactor = 3
        # Third argument is linesegment, which can be true or false. It determines
        # whether box edges around the partition elements are drawn.
        typeof(RecipesBase.apply_recipe(d, (invm, 3, true))) == Array{RecipesBase.RecipeData,1}
        typeof(RecipesBase.apply_recipe(d, (invm, 3, false))) == Array{RecipesBase.RecipeData,1}
    end

end
