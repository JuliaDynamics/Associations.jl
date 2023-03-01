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

function estimate(measure::CMIShannon, est::MesnerShalisi, x, y, z)
    e = measure.e
    (; k, w, metric) = est
    X = StateSpaceSet(x)
    Y = StateSpaceSet(y)
    Z = StateSpaceSet(z)
    joint = StateSpaceSet(X, Y, Z)
    XZ = StateSpaceSet(X, Z)
    YZ = StateSpaceSet(Y, Z)
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
        condmi += digamma(k̂)
        # Simulate ≤ by adding smallest possible nudge.
        condmi -= log(inrangecount(tree_xz, XZ[i], dmax + eps()))
        condmi -= log(inrangecount(tree_yz, YZ[i], dmax + eps()))
        condmi += log(inrangecount(tree_z, Z[i], dmax + eps()))
    end
    # The "unit" is nats.
    condmi /= N
    return _convert_logunit(condmi, ℯ, base)
end
