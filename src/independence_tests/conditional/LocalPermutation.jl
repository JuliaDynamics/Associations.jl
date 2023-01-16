using Random: shuffle!, MersenneTwister
import Statistics: quantile

export LocalPermutation
export LocalPermutationTest
export pvalue

"""
    LocalPermutation <: ConditionalIndependenceTest
    LocalPermutation(;
        measure = CMIShannon(),
        est = FrenzelPompeVelmejkaPalus(k = 5),
        kperm::Int = 5,
        nsurr::Int = 100,
        rng = Random.MersenneTwister(1234),
    )

A generic implementation of the `LocalPermutation` test  (Runge, 2018)[^Runge2018],
which tests whether two variables `X` and `Y` are conditionally independendent given a
third variable `Z` (all of which may be multivariate).

Our implementation is completely generic, so you can use any valid combination of
`definition`, `measure` and `est` that yields a conditional dependence measure
``\\hat{M}(X; Y | Z)`` with the same ordering as [`condmutualinfo`](@ref) (i.e. the conditional
variable (`Z`) is the third input variable). If a certain measure doesn't have a
`definition`, it is simply ignored.

The default is to use Shannon-type conditional mutual information, ``I(X; Y | Z)``
([`CMIShannon`](@ref) with the [`FrenzelPompeVelmejkaPalus`](@ref]) estimator.

## Description

`LocalPermutation` creates permuted `X` values using a local permutation scheme that is
based on `kperm`-th nearest neighbor searches. It attempts ``(x_i^*, y_i, z_i)_{i=1}^N``
where the goal is that ``x_i^*`` are drawn without replacement, and ``x_i`` is replaced
by ``x_j`` only if ``z_i \\approx z_j``. Then, for each permuted version of `X`, it
computes the original statistic on the permuted data, keeping `Y` and `Z` fixed, i.e.
``\\hat{M}(\\hat{X}; Y | Z)``.

[^Runge2018]: Runge, J. (2018, March). Conditional independence testing based on a
    nearest-neighbor estimator of conditional mutual information. In International
    Conference on Artificial Intelligence and Statistics (pp. 938-947). PMLR.
"""
Base.@kwdef struct LocalPermutation{M, EST, R} <: ConditionalIndependenceTest
    measure::M = CMIShannon(; base = 2)
    est::EST = FrenzelPompeVelmejkaPalus(k = 5)
    rng::R = MersenneTwister(1234)
    kperm::Int = 10
    nsurr::Int = 100
end

Base.show(io::IO, test::LocalPermutation) = print(io,
    """
    `LocalPermutation` independence test.
    -------------------------------------
    measure:        $(test.measure)
    estimator:      $(test.est)
    rng:            $(test.rng)
    # permutations: $(test.nsurr)
    k (perm)        $(test.kperm)
    """
)

"""
    LocalPermutationTest(M, Msurr, pvalue)

Holds the result of a [`LocalPermutationTest`](@ref). `M` is the measure computed on
the original data. `Msurr` is a vector of the measure computed on permuted data, where
Msurr[i] corresponds to the `i`-th permutation. `pvalue` is the `p`-value for the test.
"""
struct LocalPermutationTest{M, MS, P}
    M::M
    Msurr::MS
    pvalue::P
    nsurr::Int
end
pvalue(r::LocalPermutationTest) = r.pvalue
quantile(r::LocalPermutationTest, q) = quantile(r.Msurr, q)

function Base.show(io::IO, test::LocalPermutationTest)
    α005 = pvalue(test) < 0.05 ?
        "H₀ rejected at α = 0.05:  Yes ✅" :
        "H₀ rejected at α = 0.05:  No  ❌"
    α001 = pvalue(test) < 0.01 ?
        "H₀ rejected at α = 0.01:  Yes ✅" :
        "H₀ rejected at α = 0.01:  No  ❌"
    α0001 = pvalue(test) < 0.001 ?
        "H₀ rejected at α = 0.001: Yes ✅" :
        "H₀ rejected at α = 0.001: No  ❌"

    print(io,
        """\
        `LocalPermutation` independence test
        ----------------------------------------------------------------------------------
        H₀: "The first two variables are conditionally independent given the third"
        ----------------------------------------------------------------------------------
        Estimated: $(test.M)
        Ensemble quantiles ($(test.nsurr) permutations):
          (99.9%): $(quantile(test.Msurr, 0.999))
          (99%):   $(quantile(test.Msurr, 0.99))
          (95%):   $(quantile(test.Msurr, 0.95))
        p-value:   $(test.pvalue)
          $α005
          $α001
          $α0001\
        """

        )
end

# It is possible to specialize on the measure, e.g. LocalPermutation{CMI}. This
# should be done for the NN-based CMI methods, so we don't have to reconstruct
# KD-trees and do marginal searches for all marginals all the time.
function conditional_independence(test::LocalPermutation, x, y, z)
    (; measure, est, rng, kperm, nsurr) = test
    X, Y, Z = Dataset(x), Dataset(y), Dataset(z)
    e = test.measure.e
    @assert length(X) == length(Y) == length(Z)
    N = length(x)
    Î = estimate(measure,est, X, Y, Z)
    tree_z = KDTree(Z, Chebyshev())
    idxs_z = bulkisearch(tree_z, Z, NeighborNumber(kperm), Theiler(0))
    𝒩 = MVector{kperm, Int16}.(idxs_z) # A statically sized copy
    n̂ = collect(1:N)
    X̂ = deepcopy(X)
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
            X̂.data[i] = X.data[j]
        end
        Îs[b] = estimate( measure, est, X̂, Y, Z)
    end
    p = count(Î .<= Îs) / nsurr

    return LocalPermutationTest(Î, Îs, p, nsurr)
end

new_permutation!(n̂, rng) = shuffle!(rng, n̂)
function shuffle_neighbor_indices!(idxs::Vector{MVector{D, I}}, rng) where {D, I}
    for i = 1:length(idxs)
        shuffle!(rng, idxs[i])
    end
end
