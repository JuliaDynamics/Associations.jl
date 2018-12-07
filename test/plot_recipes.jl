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


    #@testset "StateSpaceReconstruction.AbstractEmbedding" begin
    #    @test typeof(RecipesBase.apply_recipe(d, E)) == Array{RecipesBase.RecipeData,1}
    #end

    @testset "PerronFrobenius.RectangularPartitionTransferOperator" begin
        @test typeof(RecipesBase.apply_recipe(d, invm.transfermatrix)) == Array{RecipesBase.RecipeData,1}
    end

    @testset "Invariant distribution" begin
        @test typeof(RecipesBase.apply_recipe(d, invm.measure)) == Array{RecipesBase.RecipeData,1}
    end

    @testset "RectangularInvariantMeasure" begin
        boxfillfactor = 3
        # Third argument is linesegment, which can be true or false. It determines
        # whether box edges around the partition elements are drawn.
        @test typeof(RecipesBase.apply_recipe(d, (invm, 3, true))) == Array{RecipesBase.RecipeData,1}
        @test typeof(RecipesBase.apply_recipe(d, (invm, 3, false))) == Array{RecipesBase.RecipeData,1}
    end

	@testset "Vizualizing triangulation and simplices" begin
		pts = rand(3, 15)
		n_pts = size(pts, 2)

		# Embed the point
		E = StateSpaceReconstruction.embed(pts)

		# Make sure last point is inside the convex hull of the previous points.
		# If it is not, move it towards the center of the embedding until it is.
		E_linearlyinvariant = invariantize(E)

		# Pick out points and their images (so that they can
		# be indexed with the same simplex index)
		originalpts = E.points[:, 1:end-1]
		forwardpts  = E.points[:, 2:end]

		# Re-create the embedding, but exclude the last point.
		E = StateSpaceReconstruction.embed(originalpts)

		# Triangulate all points but the last point.
		DT = delaunay(E);

		@test typeof(RecipesBase.apply_recipe(d, (pts, E, DT))) == Array{RecipesBase.RecipeData,1}
	end

end
