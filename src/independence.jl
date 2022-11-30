using Random: shuffle!, default_rng

export LocalPermutation
export independence
export ConditionalIndependenceTest

abstract type IndependenceTest end
"""
    ConditionalIndependenceTest <: IndependenceTest

The supertype for all conditional independence tests.

Currently, the following concrete implementations exist:

- [`LocalPermutation`](@ref).
"""
abstract type ConditionalIndependenceTest <: IndependenceTest end

"""
    independence(test, x, y, z)

Compute the conditional independence of `x` and `y` given `z` using
the provided [`ConditionalIndependenceTest`](@ref).
"""
function independence(test, args...; kwargs...)
    error("No concrete implementation for $(typeof(test)) test yet")
end

"""
    LocalPermutation <: ConditionalIndependenceTest
    LocalPermutation(;
        measure = CMI(),
        e::Entropy = Shannon(; base = 2),
        est = VejmelkaPalus(k = 5),
        kperm::Int = 5,
        nsurr::Int = 100,
    )

The `LocalPermutation` test  (Runge, 2018)[^Runge2018] tests whether two variables
`X` and `Y` are conditionally independendent given a third variable `Z` (all of which
may be multivariate).

The test can be used for any independence `measure`, here denoted ``\\hat{M}(X; Y | Z)``.
For example, you can use conditional mutual information, ``I(X; Y | Z)`` ([`cmi`](@ref).

## Description

`LocalPermutation` creates permuted `X` values using a local permutation scheme that is
based on `kperm`-th nearest neighbor searches. It attempts ``(x_i^*, y_i, z_i)_{i=1}^N`` where the goal is that
``x_i^*`` are drawn without replacement, and ``x_i`` is replaced by ``x_j``
only if ``z_i \\approx z_j``.

Our implementation is completely generic, so you can use any valid combination of
`method`, entropy type `e` and estimator `est` to compute the CMI.
Specialized performance-increasing methods may have been implemented for some estimators.

## Examples

```julia
using CausalityTools
using Random
x, y, z = Dataset(rand(1000, 2)), Dataset(rand(1000, 1)), Dataset(rand(1000, 1))
e = Shannon(; base = ℯ)
test_vf = LocalPermutation(; measure = CMI(MI2()), e, est = VisitationFrequency(5))
test_kr = LocalPermutation(; measure = CMI(MI2()), e, est = Kraskov(k = 10))
test_ksg1 = LocalPermutation(; measure = CMI(MI2()), e, est = KSG1(k = 10))
test_vp = LocalPermutation(; measure = CMI(), e, est = VejmelkaPalus(k = 10))

pval_vf, Îxyz_vf = independence(test_vf, x, y, z)
pval_kr, Îxyz_kr = independence(test_kr, x, y, z)
pval_ksg1, Îxyz_ksg1 = independence(test_ksg1, x, y, z)
pval_vp, Îxyz_vp = independence(test_vp, x, y, z)
```

[^Runge2018]: Runge, J. (2018, March). Conditional independence testing based on a
    nearest-neighbor estimator of conditional mutual information. In International
    Conference on Artificial Intelligence and Statistics (pp. 938-947). PMLR.
"""
Base.@kwdef struct LocalPermutation{M <: CMI, E <: Entropy, EST, R} <: ConditionalIndependenceTest
    measure::M = CMI()
    e::E = Shannon(; base = 2)
    est::EST = VejmelkaPalus(k = 5)
    rng::R = default_rng()
    kperm::Int = 5
    nsurr::Int = 200
end

# It is possible to specialize on the measure, e.g. LocalPermutation{CMI}. This
# should be done for the NN-based CMI methods, so we don't have to reconstruct
# KD-trees and do marginal searches for all marginals all the time.
function independence(test::LocalPermutation, x, y, z) # TODO: this is probably only valid for info measures.
    (; measure, e, est, rng, kperm, nsurr) = test
    @assert length(x) == length(y) == length(z)
    N = length(x)

    Î = estimate(measure, e, est, x, y, z)
    e.q ≈ 1 || throw(ArgumentError(
        "Renyi entropy with q = $(e.q) not implemented for $(typeof(est)) estimators"
    ))

    Z = Dataset(z)
    tree_z = KDTree(Z, Chebyshev())
    idxs_z = bulkisearch(tree_z, Z, NeighborNumber(kperm), Theiler(0))
    𝒩 = MVector{kperm, Int16}.(idxs_z) # A statically sized copy
    n̂ = collect(1:N)
    x̂ = deepcopy(x)
    𝒰 = zeros(Int, N) # used indices

    Îs = zeros(nsurr)
    for b in 1:nsurr
        shuffle_neighbor_indices!(𝒩, rng)
        # By re-filling, we avoid allocating extra vector for each surr. By filling with
        # zeros, we make sure that the while loop below isn't affected.
        𝒰 .= 0
        Π = new_permutation!(n̂, rng)
        for i in Π # for every point xᵢ.
            𝒩ᵢ = 𝒩[i] # shuffled neighbors to xᵢ, in terms of z
            j = first(𝒩ᵢ)
            m = 1
            while j ∈ 𝒰 && m < kperm
                m += 1
                j = 𝒩ᵢ[m]
            end
            𝒰[i] = j
            push!(𝒰, j)
            x̂[i] = x.data[j]
        end
        Îs[b] = estimate(measure, e, est, x̂, y, z)
    end
    p = count(Îs .>= Î) / nsurr
    return Î, p
end

new_permutation!(n̂, rng) = shuffle!(rng, n̂)
function shuffle_neighbor_indices!(idxs::Vector{MVector{D, I}}, rng) where {D, I}
    for i = 1:length(idxs)
        shuffle!(rng, idxs[i])
    end
end

"""
    GlobalPermutation <: ConditionalIndependenceTest

The `GlobalPermutation` conditional independence test uses surrogate
data to shuffle.

A related test is [`LocalPermutation`](@ref), but for that test, shuffled
data preserve the local

See also:
[TimeseriesSurrogates.jl](https://github.com/JuliaDynamics/TimeseriesSurrogates.jl)
"""
Base.@kwdef struct GlobalPermutation{M <: CMI, E <: Entropy, EST, R, S} <: ConditionalIndependenceTest
    measure::M = CMI()
    e::E = Shannon(; base = 2)
    est::EST = VejmelkaPalus(k = 5)
    rng::R = Random.default_rng()
    surrogate::S = 5
    nsurr::Int = 100
end
