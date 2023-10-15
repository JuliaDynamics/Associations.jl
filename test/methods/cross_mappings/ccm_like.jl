using Test
using CausalityTools
using StateSpaceSets: StateSpaceSet
n = 1000
x, y, z, w = rand(n), rand(n), StateSpaceSet(rand(n, 3)), StateSpaceSet(rand(n + 1, 3))
τ = -1

# V2.x
@testset "ConvergentCrossMapping" begin
    def = CCM(; τ)
    @test ConvergentCrossMapping() isa ConvergentCrossMapping
    @test ConvergentCrossMapping() isa CrossmapMeasure
    @test CCM() isa ConvergentCrossMapping
    @test crossmap(ExpandingSegment(def, libsizes = 100), x, y) isa Real
    @test crossmap(RandomSegment(def, libsizes = 100), x, y) isa Real
    @test crossmap(RandomVectors(def, libsizes = 100, replace = false), x, y) isa Real
    @test crossmap(RandomVectors(def, libsizes = 100, replace = true), x, y) isa Real
    @test crossmap(ExpandingSegment(def, libsizes = 100:100:500), x, y) isa Vector{<:Real}
    @test crossmap(RandomSegment(def, libsizes = 100:100:500), x, y) isa Vector{<:Real}
    @test crossmap(RandomVectors(def, libsizes = 100:100:500, replace = false), x, y) isa Vector{<:Real}
    @test crossmap(RandomVectors(def, libsizes = 100:100:500, replace = true), x, y) isa Vector{<:Real}
    @test_throws ArgumentError crossmap(RandomSegment(def, libsizes = 100), x, w) isa Real
    @test_throws ArgumentError crossmap(RandomVectors(def, libsizes = 100), x, w) isa Real

    # Ensemble analysis
    libsizes = 50
    e = Ensemble(RandomVectors(def; libsizes), nreps = 7)
    @test crossmap(e, x, y) isa Vector{<:Real}
    @test crossmap(e, x, y) |> length == 7

    libsizes = 20:10:40
    e = Ensemble(RandomVectors(def; libsizes), nreps = 7)
    @test crossmap(e, x, y) isa Vector{Vector{T}} where T
    @test crossmap(e, x, y) |> length == length(libsizes)
    @test all(length.(crossmap(e, x, y)) .== 7)
end

@testset "PairwiseAsymmetricInference" begin
    @test PairwiseAsymmetricInference() isa PairwiseAsymmetricInference
    @test PairwiseAsymmetricInference() isa CrossmapMeasure
    @test PAI() isa PairwiseAsymmetricInference

    def = PAI(; τ)

    @test crossmap(ExpandingSegment(def; libsizes = 100), x, y) isa Real
    @test crossmap(RandomSegment(def; libsizes = 100), x, y) isa Real
    @test crossmap(RandomVectors(def; libsizes = 100, replace = false), x, y) isa Real
    @test crossmap(RandomVectors(def; libsizes = 100, replace = true), x, y) isa Real
    @test crossmap(ExpandingSegment(def; libsizes = 100:100:500), x, y) isa Vector{<:Real}
    @test crossmap(RandomSegment(def; libsizes = 100:100:500), x, y) isa Vector{<:Real}
    @test crossmap(RandomVectors(def; libsizes = 100:100:500, replace = false), x, y) isa Vector{<:Real}
    @test crossmap(RandomVectors(def; libsizes = 100:100:500, replace = true), x, y) isa Vector{<:Real}
    @test_throws ArgumentError crossmap(RandomSegment(def; libsizes = 100), x, w) isa Real
    @test_throws ArgumentError crossmap(RandomVectors(def; libsizes = 100), x, w) isa Real
end

@testset "Estimator specifics" begin
    x, y = rand(50), rand(50)
    def = CCM(; τ)
    @test ExpandingSegment(def; libsizes = 100) isa CrossmapEstimator
    @test RandomSegment(def; libsizes = 100) isa CrossmapEstimator
    @test RandomVectors(def; libsizes = 100) isa CrossmapEstimator
    @test Ensemble(RandomVectors(def; libsizes = 100)) isa Ensemble
    @test crossmap(RandomVectors(def; libsizes = 10), x, y) |> length == 1
    @test crossmap(RandomVectors(def; libsizes = 10:10:20), x, y) |> length == 2
    @test crossmap(RandomVectors(def; libsizes = [10, 20, 30]), x, y) |> length == 3
    @test crossmap(RandomSegment(def; libsizes = 10), x, y) |> length == 1
    @test crossmap(RandomSegment(def; libsizes = 10:10:20), x, y) |> length == 2
    @test crossmap(RandomSegment(def; libsizes = [10, 20, 30]), x, y) |> length == 3
    @test crossmap(ExpandingSegment(def; libsizes = 10), x, y) |> length == 1
    @test crossmap(ExpandingSegment(def; libsizes = 10:10:20), x, y) |> length == 2
    @test crossmap(ExpandingSegment(def; libsizes = [10, 20, 30]), x, y) |> length == 3
end
