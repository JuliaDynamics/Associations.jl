@testset "Convergent cross mapping" begin
    @testset "Uncertainty measures" begin
        x, y = rand(100), rand(100)
        L = 50:10:100
        @test convergentcrossmap(x, y, L, uncertainty_measure = :quantile) isa Tuple{Vector{Float64}, Array{Float64, 2}}
        @test convergentcrossmap(x, y, L, uncertainty_measure = :std) isa Tuple{Vector{Float64}, Vector{Float64}}
        @test convergentcrossmap(x, y, L, quantiles = [0.25, 0.40, 0.75]) isa Tuple{Vector{Float64}, Array{Float64, 2}}
        @test convergentcrossmap(x, y, L, uncertainty_measure = :none) isa Vector{Float64}
    end

    @testset "Average measures" begin
        x, y = rand(100), rand(100)
        L = 50:10:100
        @test convergentcrossmap(x, y, L, average_measure = :median) isa Tuple{Vector{Float64}, Array{Float64, 2}}
        @test convergentcrossmap(x, y, L, average_measure = :mean) isa Tuple{Vector{Float64}, Array{Float64, 2}}
        @test convergentcrossmap(x, y, L, average_measure = :none) isa Array{Float64, 2}
    end

    @testset "Summarise" begin
        x, y = rand(100), rand(100)
        L = 50:10:100
        @test convergentcrossmap(x, y, L, summarise = true) isa Tuple{Vector{Float64}, Array{Float64, 2}}
        no_summary = convergentcrossmap(x, y, L, summarise = false)
        summary_one = no_summary[1]
        # If more than 1% of the elements in the summary of any time series length is
        # exactly 0, then assume that summary values have not been assigned.
        @test count(i->(i == 0), summary_one) < 0.99*length(summary_one)
        @test no_summary isa Vector{Vector{Float64}}
    end

    @testset "Validation functions" begin
        uncertainty_measure = :nonsense
        average_measure = :nonsense
        summarise = true
        @test_throws DomainError validate_uncertainty_measure(uncertainty_measure)
        @test_throws DomainError validate_average_measure(average_measure)
        @test_throws ErrorException validate_output_selection(:none, :none, summarise)
    end
end
