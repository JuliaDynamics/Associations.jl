
using Neighborhood: bulkisearch, inrangecount
using Neighborhood: Theiler, NeighborNumber, KDTree, Chebyshev
using SpecialFunctions: digamma

export FPVP

"""
    FPVP <: ConditionalMutualInformationEstimator
    FPVP(k = 1, w = 0)

The Frenzel-Pompe-Vejmelka-Paluš (or `FPVP` for short) estimator is used to estimate the
differential conditional mutual information using a `k`-th nearest neighbor approach that is
analogous to that of the [`KraskovStögbauerGrassberger1`](@ref) mutual information estimator
(Frenzel & Pompe, 2007[^Frenzel2007]; Vejmelka & Paluš, 2008[^Vejmelka2008]).

`w` is the Theiler window, which controls the number of temporal neighbors that are excluded
during neighbor searches.

[^Frenzel2007]:
    Frenzel, S., & Pompe, B. (2007). Partial mutual information for coupling analysis of
    multivariate time series. Physical review letters, 99(20), 204101.
    `w` is the Theiler window.
[^Vejmelka2008]:
    Vejmelka, M., & Paluš, M. (2008). Inferring the directionality of coupling with
    conditional mutual information. Physical Review E, 77(2), 026214.
"""
Base.@kwdef struct FPVP{MJ, MM} <: ConditionalMutualInformationEstimator
    k::Int = 1
    w::Int = 0
    metric_joint::MJ = Chebyshev()
    metric_marginals::MM = Chebyshev()
end

function estimate(measure::CMIShannon, est::FPVP, x, y, z)
    e = measure.e
    (; k, w, metric_joint, metric_marginals) = est
    # Ensures that vector-valued inputs are converted to StateSpaceSets, so that
    # building the marginal/joint spaces and neighbor searches are fast.
    X = StateSpaceSet(x)
    Y = StateSpaceSet(y)
    Z = StateSpaceSet(z)
    @assert length(X) == length(Y) == length(Z)
    N = length(X)
    joint = StateSpaceSet(X, Y, Z)
    XZ = StateSpaceSet(X, Z)
    YZ = StateSpaceSet(Y, Z)

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
