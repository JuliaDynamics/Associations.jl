using Test
using CausalityTools
using StateSpaceSets: Dataset
x, y, z, w = rand(1000), rand(1000), Dataset(rand(1000, 3)), Dataset(rand(1001, 3))

@testset "Measures" begin
    @test CCM() isa ConvergentCrossMapping
    @test PairwiseAsymmetricInference() isa PairwiseAsymmetricInference
end

@testset "Basic" begin
    @test crossmap(CCM(), x, y) isa Real
    @test crossmap(PAI(), x, y) isa Real
    @test_throws AssertionError crossmap(PAI(), x, w)
end

@test crossmap(CCM(), ExpandingSegment(libsizes = 100), x, y) isa Real
@test crossmap(CCM(), RandomSegment(libsizes = 100), x, y) isa Real
@test crossmap(CCM(), RandomVectors(libsizes = 100, replace = false), x, y) isa Real
@test crossmap(CCM(), RandomVectors(libsizes = 100, replace = true), x, y) isa Real
@test crossmap(CCM(), ExpandingSegment(libsizes = 100:100:500), x, y) isa Vector{<:Real}
@test crossmap(CCM(), RandomSegment(libsizes = 100:100:500), x, y) isa Vector{<:Real}
@test crossmap(CCM(), RandomVectors(libsizes = 100:100:500, replace = false), x, y) isa Vector{<:Real}
@test crossmap(CCM(), RandomVectors(libsizes = 100:100:500, replace = true), x, y) isa Vector{<:Real}
@test_throws ArgumentError crossmap(CCM(), RandomSegment(libsizes = 100), x, w) isa Real
@test_throws ArgumentError crossmap(CCM(), RandomVectors(libsizes = 100), x, w) isa Real

@test crossmap(PAI(), ExpandingSegment(libsizes = 100), x, y) isa Real
@test crossmap(PAI(), RandomSegment(libsizes = 100), x, y) isa Real
@test crossmap(PAI(), RandomVectors(libsizes = 100, replace = false), x, y) isa Real
@test crossmap(PAI(), RandomVectors(libsizes = 100, replace = true), x, y) isa Real
@test crossmap(PAI(), ExpandingSegment(libsizes = 100:100:500), x, y) isa Vector{<:Real}
@test crossmap(PAI(), RandomSegment(libsizes = 100:100:500), x, y) isa Vector{<:Real}
@test crossmap(PAI(), RandomVectors(libsizes = 100:100:500, replace = false), x, y) isa Vector{<:Real}
@test crossmap(PAI(), RandomVectors(libsizes = 100:100:500, replace = true), x, y) isa Vector{<:Real}
@test_throws ArgumentError crossmap(PAI(), RandomSegment(libsizes = 100), x, w) isa Real
@test_throws ArgumentError crossmap(PAI(), RandomVectors(libsizes = 100), x, w) isa Real

# TODO: remove for v2.0
@testset "Compat" begin
    x, y = rand(100), rand(100)
    @test crossmap(x, y, 2, 1) isa Real
    @test pai(x, y, 2, 1) isa Real
    @test length(crossmap(x, y, 2, 1, :segment, nreps = 100)) == 100
    @test length(pai(x, y, 2, 1, :segment, nreps = 100)) == 100
    @test length(crossmap(x, y, 2, 1, :random, nreps = 100)) == 100
    @test length(pai(x, y, 2, 1, :random, nreps = 100)) == 100
end
