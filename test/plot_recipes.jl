

@testset "Plotting recipes" begin
    emb = embed([diff(rand(30)) for i = 1:3], [1, 2, 3], [1, 0, -1])
    emb_invariant = invariantize(emb)

    # Test plot recipes by calling RecipesBase.apply_recipe with empty dict.
    # It should return a vector of RecipesBase.RecipeData
    d = Dict{Symbol,Any}()

    @test typeof(RecipesBase.apply_recipe(d, emb)) == Array{RecipesBase.RecipeData,1}
    @test typeof(RecipesBase.apply_recipe(d, emb_invariant)) == Array{RecipesBase.RecipeData,1}

end
