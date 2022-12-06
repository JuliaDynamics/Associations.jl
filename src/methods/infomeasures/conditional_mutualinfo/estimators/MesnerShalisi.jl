export MesnerShalisi

"""
    MesnerShalisi <: ConditionalMutualInformationEstimator
    MesnerShalisi(k = 1, w = 0)

The `MesnerShalisi` estimator is an estimator for conditional mutual information for data
that can be mixtures of discrete and continuous data (Mesner & Shalisi et al.,
2020)[^MesnerShalisi2020].

[^MesnerShalisi2020]:
    Mesner, O. C., & Shalizi, C. R. (2020). Conditional mutual information estimation for
    mixed, discrete and continuous data. IEEE Transactions on Information Theory, 67(1),
    464-484.
"""
Base.@kwdef struct MesnerShalisi{M} <: ConditionalMutualInformationEstimator
    k::Int = 1
    w::Int = 0
    metric::M = Chebyshev()
end

function estimate(infomeasure::CMI{Nothing}, e::Renyi, est::FrenzelPompe, X, Y, Z)
    e.q ≈ 1 || throw(ArgumentError(
        "Renyi entropy with q = $(e.q) not implemented for $(typeof(est)) estimators"
    ))
    (; k, w, metric) = est
    joint = Dataset(X, Y, Z)
    XZ = Dataset(X, Z)
    YZ = Dataset(Y, Z)
    Z = Dataset(Z)
    N = length(joint)
    M = 3
    tree_joint = KDTree(joint, metric)
    ds_joint = last.(bulksearch(tree_joint, joint, NeighborNumber(k), Theiler(w))[2])
    tree_xz = KDTree(XZ, metric)
    tree_yz = KDTree(YZ, metric)
    tree_z = KDTree(Z, metric)

    condmi = 0.0
    for i = 1:N
        # The notation for ρ_{i, xy} in the paper in unclear. They claim in the paper that
        # the estimator reduces to the KSG1 estimator when k̂ == k. Therefore,
        # I assume ρ_{i, xy} is the distance in the *joint* space.
        # TODO: this might not be correct..."
        dmax = ds_joint[i]
        k̂ = dmax == 0 ? inrangecount(tree_joint, joint[i], 0.0) - 1  : k
        h += digamma(k̂)
        # Simulate ≤ by adding smallest possible nudge.
        h -= log(inrangecount(tree_xz, XZ[i], dmax + eps()))
        h -= log(inrangecount(tree_yz, YZ[i], dmax + eps()))
        h += log(inrangecount(tree_z, Z[i], dmax + eps()))
    end
    h /= N

    return h / log(e.base, ℯ)
end
