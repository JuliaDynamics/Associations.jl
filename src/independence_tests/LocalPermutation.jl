using Random: shuffle!, MersenneTwister
import Statistics: quantile

export LocalPermutation
export LocalPermutationTest
export pvalue

"""
    LocalPermutation <: ConditionalIndependenceTest
    LocalPermutation(;
        measure = CMI(),
        e::Entropy = Shannon(; base = 2),
        est = VejmelkaPalus(k = 5),
        kperm::Int = 5,
        nsurr::Int = 100,
        rng = Random.default_rng(),
    )

The `LocalPermutation` test  (Runge, 2018)[^Runge2018] tests whether two variables
`X` and `Y` are conditionally independendent given a third variable `Z` (all of which
may be multivariate).

The test can be used for any independence `measure`, here denoted ``\\hat{M}(X; Y | Z)``.
For example, you can use conditional mutual information, ``I(X; Y | Z)`` ([`cmi`](@ref).

## Description

`LocalPermutation` creates permuted `X` values using a local permutation scheme that is
based on `kperm`-th nearest neighbor searches. It attempts ``(x_i^*, y_i, z_i)_{i=1}^N``
where the goal is that ``x_i^*`` are drawn without replacement, and ``x_i`` is replaced
by ``x_j`` only if ``z_i \\approx z_j``.

Our implementation is completely generic, so you can use any valid combination of
`method`, entropy type `e` and estimator `est` to compute the CMI.
Specialized performance-increasing methods may have been implemented for some estimators.

## Examples

```julia
using CausalityTools
using Random

# Create a scenario where X and Y are connected through Z. Then I(X; Y | Z) = 0.
X = randn(1000)
Y = X .+ randn(1000)
Z = randn(1000) .+ 0.5*Y
x, y, z = Dataset.([X, Y, Z])

e = Shannon(; base = ℯ)
test_vf = LocalPermutation(measure = CMIShannon(base = 2), est = VisitationFrequency(5))
test_kr = LocalPermutation(measure = CMIShannon(base = 2), est = Kraskov(k = 10))
test_fpvp = LocalPermutation(measure = CMIShannon(base = 2), est = FrenzelPompeVelmejkaPalus(k = 10))

r_vf = independence(test_vf, x, y, z)
r_kr = independence(test_kr, x, y, z)
r_fpvp = independence(test_fpvp, x, y, z)
```

[^Runge2018]: Runge, J. (2018, March). Conditional independence testing based on a
    nearest-neighbor estimator of conditional mutual information. In International
    Conference on Artificial Intelligence and Statistics (pp. 938-947). PMLR.
"""
Base.@kwdef struct LocalPermutation{M <: ConditionalMutualInformation, EST, R} <: ConditionalIndependenceTest
    measure::M = CMIShannon(; base = 2)
    est::EST = FrenzelPompeVelmejkaPalus(k = 5)
    rng::R = MersenneTwister(12345678)
    kperm::Int = 5
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
    onesided_uq = [quantile(test.Msurr, q) for q in [0.95, 0.99]]

    α005 = pvalue(test) < 0.05 ?
        "H₀ not rejected at α = 0.05 → Interpretation: X ⫫ Y | Z"  :
        "H₀ rejected at α = 0.05  → Interpretation: X and Y are dependent given Z"
    α001 = pvalue(test) < 0.01 ?
        "H₀ not rejected at α = 0.01 → Interpretation: X ⫫ Y | Z"  :
        "H₀ rejected at α = 0.01  → Interpretation: X and Y are dependent given Z"
    α0001 = pvalue(test) < 0.001 ?
        "H₀ not rejected at α = 0.01 → Interpretation: X ⫫ Y | Z"  :
        "H₀ rejected at α = 0.01  → Interpretation: X and Y are dependent given Z"

    print(io,
        """
        `LocalPermutation` independence test (H₀: X ⫫ Y | Z)"
        -------------------------------------
        Estimated: $(test.M)
        # permutations: $(test.nsurr)
        Permutation ensemble quantiles:
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
function independence(test::LocalPermutation, x, y, z)
    (; measure, est, rng, kperm, nsurr) = test
    X, Y, Z = Dataset(x), Dataset(y), Dataset(z)
    e = test.measure.e
    @assert length(X) == length(Y) == length(Z)
    N = length(x)
    Î = estimate(measure,est, X, Y, Z)

    e.q ≈ 1 || throw(ArgumentError(
        "Renyi entropy with q = $(e.q) not implemented for $(typeof(est)) estimators"
    ))
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
        Îs[b] = estimate(measure, est, X̂, Y, Z)
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
