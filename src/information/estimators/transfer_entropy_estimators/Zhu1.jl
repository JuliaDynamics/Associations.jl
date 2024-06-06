using Neighborhood: KDTree, NeighborNumber, WithinRange, Theiler, Chebyshev
using Neighborhood: bulksearch, isearch
using StateSpaceSets: AbstractStateSpaceSet, StateSpaceSet
using StateSpaceSets: dimension
using SpecialFunctions: digamma
using ComplexityMeasures: maxdists, volume_minimal_rect

export Zhu1

"""
    Zhu1 <: TransferEntropyEstimator
    Zhu1(k = 1, w = 0, base = MathConstants.e)

The `Zhu1` transfer entropy estimator [Zhu2015](@cite).

Assumes that the input data have been normalized as described in [Zhu2015](@citet).
The estimator can be used both for pairwise and conditional transfer entropy.

## Description

This estimator approximates probabilities within hyperrectangles
surrounding each point `xᵢ ∈ x` using using `k` nearest neighbor searches. However,
it also considers the number of neighbors falling on the borders of these hyperrectangles.
This estimator is an extension to the entropy estimator in [Singh2003](@citet).

`w` is the Theiler window, which determines if temporal neighbors are excluded
during neighbor searches (defaults to `0`, meaning that only the point itself is excluded
when searching for neighbours).


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


## Usage

- [`information`](@ref)`(est::Zhu1, x, y, z)`.

## Example 

```julia
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 10000)
y = rand(rng, 10000) .+ x
z = rand(rng, 10000) .+ y
est = Zhu1(TEShannon(), k = 10)
information(est, x, z, y) # should be near 0 (and can be negative)
```

## Compatible definitions

- [`TEShannon`](@ref)
"""
struct Zhu1{M} <: TransferEntropyEstimator{M}
    definition::M
    k::Int
    w::Int
end

function Zhu1(definition::M = TEShannon(); k::Int = 2, w::Int = 0) where M
    k >= 2 || throw(DomainError("The number of neighbors k must be >= 2."))
    return Zhu1(definition, k, w)
end

function information(est::Zhu1{<:TEShannon}, x::AbstractVector...)
    # The Zhu1 estimator needs to keep track of the dimension of the individual
    # terms that goes into the implicit CMI computation. We could have just used
    # `h4_marginals` here, but then we wouldn't get the dimensions out of the box.
    embedding = est.definition.embedding
    S, T, T⁺, C = individual_marginals_te(embedding, x...)
    return estimate(est, S, T, T⁺, C)
end

function estimate(est::Zhu1, 
        S::AbstractStateSpaceSet, 
        T::AbstractStateSpaceSet, 
        T⁺::AbstractStateSpaceSet, 
        C::AbstractStateSpaceSet
    )
    (; definition, k, w) = est

    joint = StateSpaceSet(S, T, T⁺, C)
    ST = StateSpaceSet(S, T, C)
    TT⁺ = StateSpaceSet(T, T⁺, C)
    T = StateSpaceSet(T, C)
    DS = dimension(S)
    DT = dimension(T)
    DT⁺ = dimension(T⁺)

    # Find distances in the joint space. Then compute, for each `xᵢ ∈ joint`, the volume of
    # the minimal rectangle containing its `k` nearest neighbors (with `k` fixed).
    W = Theiler(w)
    tree_joint = KDTree(joint, Chebyshev())
    nns_joint, ds_joint = bulksearch(tree_joint, joint, NeighborNumber(k), W)
    N = length(joint)
    ds = [ds_joint[i][k] for i = 1:N]
    vJ = volumes(joint, nns_joint, N)

    # For each `xᵢ ∈ M`, where `M` is one of the marginal spaces, count the number of
    # points within distance `ds[i]` from the point. Then count, for each point in each
    # of the marginals, how many neighbors each `xᵢ` has given `ds[i]`.
    tree_ST = KDTree(ST, Chebyshev())
    tree_TT⁺ = KDTree(TT⁺, Chebyshev())
    tree_T = KDTree(T, Chebyshev())
    nns_ST    = [isearch(tree_ST, pᵢ, WithinRange(ds[i])) for (i, pᵢ) in enumerate(ST)]
    nns_TT⁺   = [isearch(tree_TT⁺, pᵢ, WithinRange(ds[i])) for (i, pᵢ) in enumerate(TT⁺)]
    nns_T     = [isearch(tree_T, pᵢ, WithinRange(ds[i])) for (i, pᵢ) in enumerate(T)]
    kST = length.(nns_ST)
    kTT⁺ = length.(nns_TT⁺)
    kT = length.(nns_T)

    # For each `xᵢ ∈ ST`, compute the volume of the minimal rectangle of its neighbors
    # falling within distance `ds[i]` of `xᵢ` (each `xᵢ` may have a different number
    # of neighbors, since we're now using absolute *distance* to find neighbors.
    vST = volumes(ST, nns_ST, N)
    vTT⁺ = volumes(TT⁺, nns_TT⁺, N)
    vT = volumes(T, nns_T, N)

    # Compute transfer entropy
    te = mean_volumes(vJ, vST, vTT⁺, vT, N) +
        mean_digamma(kST, kTT⁺, kT, k, N, DS, DT, DT⁺)
    # Convert to target unit *after* computations, which all use natural logs.
    return _convert_logunit(te, ℯ, definition.base)
end

function volumes(x::AbstractStateSpaceSet, nn_idxs, N::Int)
    T = eltype(0.0)
    volumes = zeros(T, N)
    for (i, (xᵢ, nn_idxsᵢ)) in enumerate(zip(x, nn_idxs))
        nnsᵢ = @views x[nn_idxsᵢ] # the actual coordinates of the points
        distsᵢ = maxdists(xᵢ, nnsᵢ)
        volumes[i] = volume_minimal_rect(distsᵢ)
    end
    return volumes
end

function mean_volumes(vols_joint, vols_ST, vols_TT⁺, vols_T, N::Int)
    vol = 0.0
    n_ignore = 0
    for i = 1:N
        num = vols_TT⁺[i] * vols_ST[i]
        den = vols_joint[i] * vols_T[i]
        if den != 0
            vol += log(num / den)
        else
            n_ignore += 1
        end
    end
    return vol / (N - n_ignore)
end

function mean_digamma(ks_ST, ks_TT⁺, ks_T, k::Int, N::Int,
        dS::Int, dT::Int, dT⁺::Int)

    Φ = 0.0
    for i = 1:N
        Φ += digamma(k) +
            digamma(ks_T[i]) -
            digamma(ks_TT⁺[i]) -
            digamma(ks_ST[i]) +
            (dT⁺ + dT - 1) / ks_TT⁺[i] +
            (dS + dT - 1) / ks_ST[i] -
            (dT⁺ + dT + dS - 1) / k -
            (dT - 1) / ks_T[i]
    end
    return Φ / N
end
