using Test
using CausalityTools
using StateSpaceSets: Dataset
n = 1000
x, y, z, w = rand(n), rand(n), Dataset(rand(n, 3)), Dataset(rand(n + 1, 3))

τ = -1

# Deprecated
@test crossmap(x, y, 3, τ) isa Float64
@test crossmap(x, y, 3, τ , :random) isa Vector{Float64}
@test crossmap(x, y, 3, τ, :segment) isa Vector{Float64}

# V2.x
@testset "ConvergentCrossMapping" begin
    @test ConvergentCrossMapping() isa ConvergentCrossMapping
    @test ConvergentCrossMapping() isa CrossmapMeasure
    @test CCM() isa ConvergentCrossMapping
    @test crossmap(CCM(; τ), ExpandingSegment(libsizes = 100), x, y) isa Real
    @test crossmap(CCM(; τ), RandomSegment(libsizes = 100), x, y) isa Real
    @test crossmap(CCM(; τ), RandomVectors(libsizes = 100, replace = false), x, y) isa Real
    @test crossmap(CCM(; τ), RandomVectors(libsizes = 100, replace = true), x, y) isa Real
    @test crossmap(CCM(; τ), ExpandingSegment(libsizes = 100:100:500), x, y) isa Vector{<:Real}
    @test crossmap(CCM(; τ), RandomSegment(libsizes = 100:100:500), x, y) isa Vector{<:Real}
    @test crossmap(CCM(; τ), RandomVectors(libsizes = 100:100:500, replace = false), x, y) isa Vector{<:Real}
    @test crossmap(CCM(; τ), RandomVectors(libsizes = 100:100:500, replace = true), x, y) isa Vector{<:Real}
    @test_throws ArgumentError crossmap(CCM(; τ), RandomSegment(libsizes = 100), x, w) isa Real
    @test_throws ArgumentError crossmap(CCM(; τ), RandomVectors(libsizes = 100), x, w) isa Real

    # Ensemble analysis
    libsizes = 50
    e = Ensemble(CCM(), RandomVectors(; libsizes), nreps = 7)
    @test crossmap(e, x, y) isa Vector{<:Real}
    @test crossmap(e, x, y) |> length == 7

    libsizes = 20:10:40
    e = Ensemble(CCM(), RandomVectors(; libsizes), nreps = 7)
    @test crossmap(e, x, y) isa Vector{Vector{T}} where T
    @test crossmap(e, x, y) |> length == length(libsizes)
    @test all(length.(crossmap(e, x, y)) .== 7)
end

@testset "PairwiseAsymmetricInference" begin
    @test PairwiseAsymmetricInference() isa PairwiseAsymmetricInference
    @test PairwiseAsymmetricInference() isa CrossmapMeasure
    @test PAI() isa PairwiseAsymmetricInference

    
    @test crossmap(PAI(; τ), ExpandingSegment(libsizes = 100), x, y) isa Real
    @test crossmap(PAI(; τ), RandomSegment(libsizes = 100), x, y) isa Real
    @test crossmap(PAI(; τ), RandomVectors(libsizes = 100, replace = false), x, y) isa Real
    @test crossmap(PAI(; τ), RandomVectors(libsizes = 100, replace = true), x, y) isa Real
    @test crossmap(PAI(; τ), ExpandingSegment(libsizes = 100:100:500), x, y) isa Vector{<:Real}
    @test crossmap(PAI(; τ), RandomSegment(libsizes = 100:100:500), x, y) isa Vector{<:Real}
    @test crossmap(PAI(; τ), RandomVectors(libsizes = 100:100:500, replace = false), x, y) isa Vector{<:Real}
    @test crossmap(PAI(; τ), RandomVectors(libsizes = 100:100:500, replace = true), x, y) isa Vector{<:Real}
    @test_throws ArgumentError crossmap(PAI(; τ), RandomSegment(libsizes = 100), x, w) isa Real
    @test_throws ArgumentError crossmap(PAI(; τ), RandomVectors(libsizes = 100), x, w) isa Real
end

@testset "Estimator specifics" begin
    x, y = rand(50), rand(50)
    @test ExpandingSegment(libsizes = 100) isa CrossmapEstimator
    @test RandomSegment(libsizes = 100) isa CrossmapEstimator
    @test RandomVectors(libsizes = 100) isa CrossmapEstimator
    @test Ensemble(CCM(), RandomVectors(libsizes = 100)) isa Ensemble
    @test crossmap(PAI(; τ), RandomVectors(libsizes = 10), x, y) |> length == 1
    @test crossmap(PAI(; τ), RandomVectors(libsizes = 10:10:20), x, y) |> length == 2
    @test crossmap(PAI(; τ), RandomVectors(libsizes = [10, 20, 30]), x, y) |> length == 3
    @test crossmap(PAI(; τ), RandomSegment(libsizes = 10), x, y) |> length == 1
    @test crossmap(PAI(; τ), RandomSegment(libsizes = 10:10:20), x, y) |> length == 2
    @test crossmap(PAI(; τ), RandomSegment(libsizes = [10, 20, 30]), x, y) |> length == 3
    @test crossmap(PAI(; τ), ExpandingSegment(libsizes = 10), x, y) |> length == 1
    @test crossmap(PAI(; τ), ExpandingSegment(libsizes = 10:10:20), x, y) |> length == 2
    @test crossmap(PAI(; τ), ExpandingSegment(libsizes = [10, 20, 30]), x, y) |> length == 3
end

# TODO: remove for v2.0
@testset "Compat" begin
    x, y = rand(100), rand(100)
    @test crossmap(x, y, 2, τ) isa Real
    @test pai(x, y, 2, τ) isa Real
    @test length(crossmap(x, y, 2, τ, :segment, nreps = 100)) == 100
    @test length(pai(x, y, 2, τ, :segment, nreps = 100)) == 100
    @test length(crossmap(x, y, 2, τ, :random, nreps = 100)) == 100
    @test length(pai(x, y, 2, τ, :random, nreps = 100)) == 100
end
