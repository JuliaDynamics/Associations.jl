using SpecialFunctions: digamma
using Neighborhood: KDTree, NeighborNumber, WithinRange, Euclidean
using Neighborhood: bulksearch, isearch

export Lindner

"""
    Lindner <: TransferEntropyEstimator
    Lindner(definition = Shannon(); k = 1, w = 0, base = 2)

The `Lindner` transfer entropy estimator [Lindner2011](@cite), which is
also used in the Trentool MATLAB toolbox, and is based on nearest neighbor searches.

## Usage

- Use with [`association`](@ref) to compute [`TEShannon`](@ref) from input data.

## Keyword parameters

`w` is the Theiler window, which determines if temporal neighbors are excluded
during neighbor searches (defaults to `0`, meaning that only the point itself is excluded
when searching for neighbours).

The estimator can be used both for pairwise and conditional transfer entropy estimation.

## Description

For a given points in the joint embedding space `jᵢ`, this estimator first computes the
distance `dᵢ` from `jᵢ` to its `k`-th nearest neighbor. Then, for each point `mₖ[i]` in
the `k`-th marginal space, it counts the number of points within radius `dᵢ`.

The Shannon transfer entropy is then computed as

```math
TE_S(X \\to Y) =
\\psi(k) + \\dfrac{1}{N} \\sum_{i}^n
\\left[
    \\sum_{k=1}^3 \\left( \\psi(m_k[i] + 1) \\right)
\\right],
```

where the index `k` references the three marginal subspaces `T`, `TTf` and `ST` for which
neighbor searches are performed. Here this estimator has been modified to allow for 
conditioning too (a simple modification to [Lindner2011](@citet)'s equation 5 and 6). 

## Example 

```julia
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 10000)
y = rand(rng, 10000) .+ x
z = rand(rng, 10000) .+ y
est = Lindner(TEShannon(), k = 10)
information(est, x, z, y) # should be near 0 (and can be negative)
```

## Compatible definitions

- [`TEShannon`](@ref)
"""
struct Lindner{E} <: TransferEntropyEstimator{E}
    definition::E
    k::Int # number of neighbors in joint space.
    w::Int

    function Lindner(definition::E = TEShannon(); k::Int = 2, w::Int = 0) where {E}
        k >= 2 || throw(DomainError("The number of neighbors k must be >= 2."))
        new{E}(definition, k, w)
    end
end

function information(est::Lindner{<:TEShannon}, x::AbstractVector...)
    verify_number_of_inputs_vars(est.definition, length(x))
    S, T, T⁺, C = individual_marginals_te(est.definition.embedding, x...)
    return estimate(est, S, T, T⁺, C)
end

# This method is separate from the one above because when using `SurrogateAssociationTest`,
# `S` is repeatedly shuffled, while the other marginals are not, so we avoid
# allocating a bunch of new StateSpaceSets for every shuffle.
function estimate(est::Lindner{<:TEShannon},
        S::AbstractStateSpaceSet,
        T::AbstractStateSpaceSet,
        T⁺::AbstractStateSpaceSet,
        C::AbstractStateSpaceSet)

    # This layer ensures that the number of `StateSpaceSet`s that must be 
    # constructed is minimal when doing e.g. surrogate testing (then,
    # `S` is the only marginal changing).
    TT⁺C = StateSpaceSet(T, T⁺, C)
    TC = StateSpaceSet(T, C)
    return estimate_with_premade_embeddings(est, S, TT⁺C, TC)
end

function estimate_with_premade_embeddings(est::Lindner{<:TEShannon},
        S::AbstractStateSpaceSet,
        TT⁺C::AbstractStateSpaceSet,
        TC::AbstractStateSpaceSet)
    (; definition, k, w) = est

    joint = StateSpaceSet(S, TT⁺C)
    STC = StateSpaceSet(S, TC)
    N = length(joint)
    W = Theiler(w)
    metric =  Chebyshev()
    tree_joint = KDTree(joint, metric)
    nns_joint, ds_joint = bulksearch(tree_joint, joint, NeighborNumber(k), W)
    # For each `xᵢ ∈ M`, where `M` is one of the marginal spaces, count the number of
    # points within distance `ds[i]` from the point. Then count, for each point in each
    # of the marginals, how many neighbors each `xᵢ` has given `ds[i]`.
    ds = last.(ds_joint) # only care about distance to the k-th neighbor
    tree_STC = KDTree(STC, metric)
    tree_TT⁺C = KDTree(TT⁺C, metric)
    tree_TC = KDTree(TC, metric)
    nns_STC  = [isearch(tree_STC, pᵢ, WithinRange(ds[i])) for (i, pᵢ) in enumerate(STC)]
    nns_TT⁺C = [isearch(tree_TT⁺C, pᵢ, WithinRange(ds[i])) for (i, pᵢ) in enumerate(TT⁺C)]
    nns_TC   = [isearch(tree_TC, pᵢ, WithinRange(ds[i])) for (i, pᵢ) in enumerate(TC)]

    n_STC = length.(nns_STC)
    n_TT⁺C = length.(nns_TT⁺C)
    n_TC = length.(nns_TC)
    te = 0.0
    for i = 1:N
        te += digamma(n_TC[i] + 1) - digamma(n_TT⁺C[i] + 1) - digamma(n_STC[i])
    end
    te /= N
    # The "unit" is nats
    te += digamma(k)

    # Convert to target base *after* digamma computations, because the digamma function
    # is a function of the natural log.
    return _convert_logunit(te, ℯ, definition.base)
end
