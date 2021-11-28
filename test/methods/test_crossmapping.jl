using CausalityTools, StatsBase

@testset "Cross mapping" begin

	x, y = rand(100), rand(100)
	@test crossmap(x, y, 3, 1) isa Float64
	@test crossmap(x, y, 3, 1, :random) isa Vector{Float64}
	@test crossmap(x, y, 3, 1, :segment) isa Vector{Float64}
end

@testset "Pairwise asymmetric inference" begin
	@test pai(x, y, 3, 1) isa Float64
	@test pai(x, y, 3, 1, :random) isa Vector{Float64}
	@test pai(x, y, 3, 1, :segment) isa Vector{Float64}
end
