using Random: shuffle!
using Random
import Statistics: quantile

export LocalPermutationTest
export LocalPermutationTestResult
export pvalue

# `LocalPermutationClosenessSearch`and its subtypes is just for internal use.
# The `LocalPermutationTest` as given in Runge et al. (2018) uses neighbor searches,
# for determining closeness. This limits the usefulness of the test mostly to continuous
# data. It is possible to extend the test for discrete/mixed data, but then other
# "closeness schemes" must be applied. By dispatching on `LocalPermutationClosenessSearch`
# internally, we don't need to introduce breaking changes later when other search types
# are introduced.
"""
The supertype of all types indicating a way of determining "closeness" for the
local permutation algorithm.
"""
abstract type LocalPermutationClosenessSearch end

"""
    NeighborCloseness <: LocalPermutationClosenessSearch

Determine closeness between points for [`LocalPermutationTest`](@ref) using nearest neighbors
searches.
"""
struct NeighborCloseness <: LocalPermutationClosenessSearch end

"""
    LocalPermutationTest <: IndependenceTest
    LocalPermutationTest(measure, [est];
        kperm::Int = 5,
        nshuffles::Int = 100,
        rng = Random.default_rng())

`LocalPermutationTest` is a generic conditional independence test (Runge, 2018)[^Runge2018]
for assessing whether two variables `X` and `Y` are conditionally independendent given a
third variable `Z` (all of which may be multivariate).

Any association `measure` (with a compatible estimator `est`, if relevant) with ordering
``\\hat{M}(X; Y | Z)`` (conditional variable is the third) can be used. To obtain the
nearest-neighbor approach in Runge, 2018, use the [`CMIShannon`](@ref) measure with the
[`FPVP`](@ref) estimator.

## Description

This is a generic one-sided hypothesis test that checks whether `x` and `y`
are independent (given `z`, if provided) based on resampling from a null distribution
assumed to represent independence between the variables. The null distribution is generated
by repeatedly shuffling the input data in some way that is intended
to break any dependence between the input variables.

For each shuffle, the provided `measure` is computed (using `est`,
if relevant) while keeping `Y` and `Z` fixed, but permuting `X`,
i.e. ``\\hat{M}(\\hat{X}; Y | Z)``. Each shuffle of `X` is done conditional on `Z`, such
that `xᵢ` is replaced with `xⱼ` only if `zᵢ ≈ zⱼ`, i.e. `zᵢ` and `zⱼ` are close.
Closeness is determined by a `kperm`-th nearest neighbor search among the points in `Z`,
and permuted points are constructed as
``(x_i^*, y_i, z_i)_{i=1}^N``, where the goal is that ``x_i^*`` are drawn without
replacement, and ``x_i`` is replaced by ``x_j`` only if ``z_i \\approx z_j``.
This procedure is repeated `nshuffles` times, and a test summary is returned.

## Examples

See [quickstart examples](@ref quickstart_localpermutationtest).

[^Runge2018]: Runge, J. (2018, March). Conditional independence testing based on a
    nearest-neighbor estimator of conditional mutual information. In International
    Conference on Artificial Intelligence and Statistics (pp. 938-947). PMLR.

"""
struct LocalPermutationTest{M, EST, C, R} <: IndependenceTest
    measure::M
    est::EST
    rng::R
    kperm::Int
    nshuffles::Int
    closeness_search::C
    function LocalPermutationTest(measure::M, est::EST;
            rng::R = Random.default_rng(),
            kperm::Int = 10,
            nshuffles::Int = 100,
            closeness_search::C = NeighborCloseness()) where {M, EST, C, R}
        new{M, EST, C, R}(measure, est, rng, kperm, nshuffles, closeness_search)
    end
end

Base.show(io::IO, test::LocalPermutationTest) = print(io,
    """
    `LocalPermutationTest` independence test.
    -------------------------------------
    measure:    $(test.measure)
    estimator:  $(test.est)
    rng:        $(test.rng)
    # shuffles: $(test.nshuffles)
    k (perm)    $(test.kperm)
    """
)

"""
    LocalPermutationTestResult(M, Msurr, pvalue)

Holds the result of a [`LocalPermutationTestTestResult`](@ref). `M` is the measure computed on
the original data. `Msurr` is a vector of the measure computed on permuted data, where
Msurr[i] corresponds to the `i`-th permutation. `pvalue` is the `p`-value for the test.
"""
struct LocalPermutationTestResult{M, MS, P}
    M::M
    Msurr::MS
    pvalue::P
    nshuffles::Int
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
        Ensemble quantiles ($(test.nshuffles) permutations):
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
    (; measure, est, rng, kperm, nshuffles) = test
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
    Îs = zeros(nshuffles)
    for b in 1:nshuffles
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
    p = count(Î .<= Îs) / nshuffles

    return LocalPermutationTestResult(Î, Îs, p, nshuffles)
end

new_permutation!(n̂, rng) = shuffle!(rng, n̂)
function shuffle_neighbor_indices!(idxs::Vector{MVector{D, I}}, rng) where {D, I}
    for i = 1:length(idxs)
        shuffle!(rng, idxs[i])
    end
end


include("transferentropy.jl")
