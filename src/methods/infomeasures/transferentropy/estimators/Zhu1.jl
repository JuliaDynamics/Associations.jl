
using Neighborhood: KDTree, NeighborNumber, WithinRange, Theiler, Chebyshev
using Neighborhood: bulksearch, isearch
using DelayEmbeddings: AbstractDataset
using SpecialFunctions: digamma

export Zhu1

"""
    Zhu1 <: EntropyEstimator
    Zhu1(k = 1, w = 0, base = MathConstants.e)

The `Zhu1` transfer entropy estimator (Zhu et al., 2015)[^Zhu2015].

This estimator approximates probabilities within hyperrectangles
surrounding each point `xᵢ ∈ x` using using `k` nearest neighbor searches. However,
it also considers the number of neighbors falling on the borders of these hyperrectangles.
This estimator is an extension to the entropy estimator in Singh et al. (2003).

`w` is the Theiler window, which determines if temporal neighbors are excluded
during neighbor searches (defaults to `0`, meaning that only the point itself is excluded
when searching for neighbours).

[^Zhu2015]:
    Zhu, J., Bellanger, J. J., Shu, H., & Le Bouquin Jeannès, R. (2015). Contribution to
    transfer entropy estimation via the k-nearest-neighbors approach. Entropy, 17(6),
    4173-4201.
[^Singh2003]:
    Singh, H., Misra, N., Hnizdo, V., Fedorowicz, A., & Demchuk, E. (2003). Nearest
    neighbor estimates of entropy. American journal of mathematical and management
    sciences, 23(3-4), 301-321.
"""
Base.@kwdef struct Zhu1{B} <: EntropyEstimator
    k::Int = 1
    w::Int = 0
    base::B = MathConstants.e

    function Zhu1(k::Int, w::Int, base::B) where B
        k >= 2 || throw(DomainError("The number of neighbors k must be >= 2."))
        new{B}(k, w, base)
    end
end

function te(est::Zhu1,
    S::AbstractDataset{DS, Q},
    T::AbstractDataset{DT, Q},
    T⁺::AbstractDataset{DTT, Q}
) where {DS, DT, DTT, Q}
    (; k, w, base) = est
    joint = Dataset(S, T, T⁺)
    ST, TT⁺ = Dataset(S, T), Dataset(T, T⁺)

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
    tree_ST, tree_TT⁺, tree_T = KDTree.([ST, TT⁺, T], Ref(Chebyshev()))
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
    return mean_volumes(vJ, vST, vTT⁺, vT, N) +
        mean_digamma(kST, kTT⁺, kT, k, N, DS, DT, DTT)
end

function volumes(x::AbstractDataset, nn_idxs, N::Int)
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
    for i = 1:N
        vol += (vols_TT⁺[i] * vols_ST[i]) / (vols_joint[i] * vols_T[i])
    end
    return vol / N
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