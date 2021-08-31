using StatsBase

@testset "Cross mapping" begin

	@testset "Prediction lags" begin
	    x, y = rand(100), rand(100)
	    crossmap(x, y, η = 1)
	    crossmap(x, y, η = -1)
	    crossmap(x, y, η = 5)
	    crossmap(x, y, η = -5)
	end

	@testset "Embedding params" begin
	    x, y = rand(100), rand(100)
	    @test_throws DomainError crossmap(x, y, dim = -3)
	    @test_throws DomainError crossmap(x, y, dim = 100, τ = 2)
	    crossmap(x, y, dim = 3)
	    crossmap(x, y, dim = 10, τ = 2)
	end

	@testset "Exclusion radii" begin
	    x, y = rand(100), rand(100)
	    [crossmap(x, y, theiler_window = i) for i in rand(1:25, 10)]
	    #@test_throws DomainError crossmap(x, y, theiler_window = -1)
	end

	@testset "Replacements" begin
	    x, y = rand(100), rand(100)
	    crossmap(x, y, replace = false)
	    crossmap(x, y, replace = true)
	end

	@testset "Correspondence measures" begin
	    @testset "$i" for i in 1:100
	        x, y = rand(100), rand(100)
	        @test all(crossmap(x, y, correspondence_measure = StatsBase.rmsd) .>= 0)
	        @test all([-1 <= x <= 1 for x in crossmap(x, y, correspondence_measure = StatsBase.cor)])
	    end
	end

end
