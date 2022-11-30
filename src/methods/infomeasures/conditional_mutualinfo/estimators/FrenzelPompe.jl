
using Neighborhood: bulkisearch, inrangecount
using Neighborhood: Theiler, NeighborNumber, KDTree, Chebyshev
using SpecialFunctions: digamma

export FrenzelPompe

"""
    FrenzelPompe <: ConditionalMutualInformationEstimator
    FrenzelPompe(k = 1, w = 0)

The `FrenzelPompe` estimator is used to estimate the differential conditional
mutual information using a `k`-th nearest neighbor approach that is
analogous to that of the [`KSG1`](@ref) mutual information estimator
(Frenzel & Pompe, 2007).

This estimator is identical to the [`VejmelkaPalus`](@ref) estimator,
which appeared in a separate paper around the same time.

`w` is the Theiler window, which controls the number of temporal neighbors that are excluded
during neighbor searches.

[^Frenzel2007]:
    Frenzel, S., & Pompe, B. (2007). Partial mutual information for coupling analysis of
    multivariate time series. Physical review letters, 99(20), 204101.
"""
Base.@kwdef struct FrenzelPompe{MJ, MM} <: ConditionalMutualInformationEstimator
    k::Int = 1
    w::Int = 0
    metric_joint::MJ = Chebyshev()
    metric_marginals::MM = Chebyshev()
end

function estimate(infomeasure::CMI{Nothing}, e::Renyi, est::FrenzelPompe, X, Y, Z)
    e.q ≈ 1 || throw(ArgumentError(
        "Renyi entropy with q = $(e.q) not implemented for $(typeof(est)) estimators"
    ))
    (; k, w, metric_joint, metric_marginals) = est
    @assert length(X) == length(Y) == length(Z)
    N = length(X)
    # Ensures that vector-valued inputs are converted to datasets, so that
    # building the marginal/joint spaces and neighbor searches are fast.
    X = Dataset(X)
    Y = Dataset(Y)
    Z = Dataset(Z)
    joint = Dataset(X, Y, Z)
    XZ = Dataset(X, Z)
    YZ = Dataset(Y, Z)

    tree_joint = KDTree(joint, metric_joint)
    ds_joint = last.(bulksearch(tree_joint, joint, NeighborNumber(k), Theiler(w))[2])
    tree_xz = KDTree(XZ, metric_marginals)
    tree_yz = KDTree(YZ, metric_marginals)
    tree_z = KDTree(Z, metric_marginals)

    condmi = 0.0
    for (i, dᵢ) in enumerate(ds_joint)
        # Usually, we subtract 1 because inrangecount includes the point itself,
        # but we'll have to add it again inside the digamma, so just skip it.
        condmi += digamma(k)
        condmi -= digamma(inrangecount(tree_xz, XZ[i], dᵢ))
        condmi -= digamma(inrangecount(tree_yz, YZ[i], dᵢ))
        condmi += digamma(inrangecount(tree_z, Z[i], dᵢ))
    end
    condmi /= N

    return condmi / log(e.base, ℯ)
end

estimate(infomeasure::CMI, est::FrenzelPompe, args...; base = 2, kwargs...) =
    estimate(infomeasure, Shannon(; base), est, args...; kwargs...)
