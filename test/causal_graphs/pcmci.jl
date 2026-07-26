using Associations
using StateSpaceSets: StateSpaceSet
using Random
using Distributions: Uniform

# Test on a system with known lagged couplings and no contemporaneous coupling.
# This is a linear system, so we can use `CorrTest` for quick estimation.
function linear5(npts::Int; rng=Random.default_rng(),
    a=0.5,      # autoregressive strength
    c=0.2,      # coupling strength
    σ=0.05)      # noise strength
    # `pad` extra transient points are generated (and later discarded) so that
    # the longest lag (5) always has valid history.
    pad = 100
    N = npts + pad
    x1 = zeros(N);
    x2 = zeros(N);
    x3 = zeros(N);
    x4 = zeros(N);
    x5 = zeros(N)
    u = Uniform(max(a * 0.5, 0.1), min(a*1.5, 0.9))
    for t in 6:N
        x1[t] = rand(rng, u)*x1[t-1] + σ*randn(rng)
        x2[t] = rand(rng, u)*x2[t-1] + c*x1[t-2] + σ*randn(rng)
        x3[t] = rand(rng, u)*x3[t-1] + c*x2[t-3] + σ*randn(rng)
        x4[t] = rand(rng, u)*x4[t-1] + c*x1[t-4] + σ*randn(rng)
        x5[t] = rand(rng, u)*x5[t-1] + c*x3[t-1] + c*x4[t-5] + σ*randn(rng)
    end
    return StateSpaceSet(x1[(pad+1):end], x2[(pad+1):end], x3[(pad+1):end],
        x4[(pad+1):end], x5[(pad+1):end])
end
X = linear5(800; rng=Xoshiro(1234));

# `τmax = 7` so the longest coupling lag (5) can be detected.
alg = PCMCI(; utest=CorrTest(), ctest=CorrTest(), α=0.05, τmax=7)
parents = infer_graph(alg, X, verbose=false)

@test parents isa Vector{<:PCMCISelectedParents}

# `has_link(parents, i, j, τ)` returns `true` if the driver `xⱼ(τ)` was selected as
# a parent of the target `xᵢ(0)`, i.e. if the link `xⱼ(τ) → xᵢ(0)` was identified.
function has_link(parents, i, j, τ)
    p = parents[i]
    return any(((pj, pτ),) -> pj == j && pτ == τ, zip(parents_js(p), parents_τs(p)))
end

@testset "Autoregressive self-links xᵢ(-1) → xᵢ" begin
    @test has_link(parents, 1, 1, -1)
    @test has_link(parents, 2, 2, -1)
    @test has_link(parents, 3, 3, -1)
    @test has_link(parents, 4, 4, -1)
    @test has_link(parents, 5, 5, -1)
end

@testset "Coupling links" begin
    @test has_link(parents, 2, 1, -2) # x1(-2) → x2
    @test has_link(parents, 3, 2, -3) # x2(-3) → x3
    @test has_link(parents, 4, 1, -4) # x1(-4) → x4
    @test has_link(parents, 5, 3, -1) # x3(-1) → x5
    @test has_link(parents, 5, 4, -5) # x4(-5) → x5
end
