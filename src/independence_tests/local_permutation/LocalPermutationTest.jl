using Random: shuffle!, MersenneTwister
import Statistics: quantile

export LocalPermutationTest
export LocalPermutationTestResult
export pvalue

"""
    LocalPermutationTest <: ConditionalIndependenceTest
    LocalPermutationTest(measure = CMIShannon(), est = FPVP(k = 5);
        kperm::Int = 5,
        nsurr::Int = 100,
        rng = Random.MersenneTwister(1234))

A generic implementation of the `LocalPermutationTest` test  (Runge, 2018)[^Runge2018],
which tests whether two variables `X` and `Y` are conditionally independendent given a
third variable `Z` (all of which may be multivariate).

You can use any valid combination of `definition`, `measure` and `est` that yields a
conditional dependence measure ``\\hat{M}(X; Y | Z)`` with the same ordering as
[`condmutualinfo`](@ref) (i.e. the conditional variable (`Z`) is the third input variable).
The default measure is Shannon conditional mutual information [`CMIShannon`](@ref),
``I(X; Y | Z)`` with the nearest-neighbor based [`FPVP`](@ref]) differential CMI estimator.

## Description

`LocalPermutationTest` creates permuted `X` values using a local permutation scheme that is
based on `kperm`-th nearest neighbor searches. Permuted points are constructed as
``(x_i^*, y_i, z_i)_{i=1}^N``, where the goal is that ``x_i^*`` are drawn without
replacement, and ``x_i`` is replaced by ``x_j`` only if ``z_i \\approx z_j``.
Then, for each permuted version of `X`, it computes the original statistic on the permuted
data, keeping `Y` and `Z` fixed, i.e. ``\\hat{M}(\\hat{X}; Y | Z)``.

## Example usage

```julia
x = randn(2000)
y = randn(2000) .+ 0.7 .* x
z = sin.(randn(2000)) .* 0.5 .* y
test = LocalPermutationTest(CMIShannon(; base = 2), FPVP(k = 10))
independence(test, x, y, z)
```

[^Runge2018]: Runge, J. (2018, March). Conditional independence testing based on a
    nearest-neighbor estimator of conditional mutual information. In International
    Conference on Artificial Intelligence and Statistics (pp. 938-947). PMLR.
"""
struct LocalPermutationTest{M, EST, R} <: ConditionalIndependenceTest
    measure::M
    est::EST
    rng::R
    kperm::Int
    nsurr::Int
    function LocalPermutationTest(measure::M = CMIShannon(; base = 2), est::EST = FPVP(k = 5);
            rng::R = MersenneTwister(1234),
            kperm::Int = 10,
            nsurr::Int = 100) where {M, EST, R}
        new{M, EST, R}(measure, est, rng, kperm, nsurr)
    end
end

Base.show(io::IO, test::LocalPermutationTest) = print(io,
    """
    `LocalPermutationTest` independence test.
    -------------------------------------
    measure:        $(test.measure)
    estimator:      $(test.est)
    rng:            $(test.rng)
    # permutations: $(test.nsurr)
    k (perm)        $(test.kperm)
    """
)

"""
    LocalPermutationTestResult(M, Msurr, pvalue)

Holds the result of a [`LocalPermutationTestTestResult`](@ref). `M` is the measure computed on
the original data. `Msurr` is a vector of the measure computed on permuted data, where
Msurr[i] corresponds to the `i`-th permutation. `pvalue` is the `p`-value for the test.
"""
struct LocalPermutationTestTestResult{M, MS, P}
    M::M
    Msurr::MS
    pvalue::P
    nsurr::Int
end
pvalue(r::LocalPermutationTestResult) = r.pvalue
quantile(r::LocalPermutationTestResult, q) = quantile(r.Msurr, q)

function Base.show(io::IO, test::LocalPermutationTestResult)
    α005 = pvalue(test) < 0.05 ?
        "α = 0.05: ✓ Evidence favors dependence" :
        "α = 0.05: ✖ Independence cannot be rejected"
    α001 = pvalue(test) < 0.01 ?
        "α = 0.01: ✓ Evidence favors dependence" :
        "α = 0.01: ✖ Independence cannot be rejected"
    α0001 = pvalue(test) < 0.001 ?
        "α = 0.001: ✓ Evidence favors dependence" :
        "α = 0.001: ✖ Independence cannot be rejected"

    print(io,
        """\
        `LocalPermutationTest` independence test
        ----------------------------------------------------------------------------------
        H₀: "The first two variables are conditionally independent given the 3rd variable"
        H₁: "The first two variables are conditionally dependent given the 3rd variable"
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

function independence(test::LocalPermutationTest, x, y)
    throw(ArgumentError("`LocalPermutationTest` is a conditional independence test, and thus must be given three input variables. Only two were given."))
end

# It is possible to specialize on the measure, e.g. LocalPermutationTest{CMI}. This
# should be done for the NN-based CMI methods, so we don't have to reconstruct
# KD-trees and do marginal searches for all marginals all the time.
function independence(test::LocalPermutationTest, x, y, z)
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

    return LocalPermutationTestResult(Î, Îs, p, nsurr)
end

new_permutation!(n̂, rng) = shuffle!(rng, n̂)
function shuffle_neighbor_indices!(idxs::Vector{MVector{D, I}}, rng) where {D, I}
    for i = 1:length(idxs)
        shuffle!(rng, idxs[i])
    end
end


include("transferentropy.jl")
