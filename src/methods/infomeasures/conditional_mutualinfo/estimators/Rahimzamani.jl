export Rahimzamani

"""
    Rahimzamani <: ConditionalMutualInformationEstimator
    Rahimzamani(k = 1, w = 0)

The `Rahimzamani` estimator, short for Rahimzamani-Asnani-Viswanath-Kannan,
is an estimator for conditional mutual information for data that can be mixtures of
discrete and continuous data (Rahimzamani et al., 2018)[^Rahimzamani2018].

This is very similar to the [`GaoKannanOhViswanath`](@ref) mutual information estimator,
but has been expanded to the conditional case.


[^Rahimzamani2018]:
    Rahimzamani, A., Asnani, H., Viswanath, P., & Kannan, S. (2018). Estimators for
    multivariate information measures in general probability spaces. Advances in Neural
    Information Processing Systems, 31.
"""
Base.@kwdef struct Rahimzamani{M} <: ConditionalMutualInformationEstimator
    k::Int = 1
    w::Int = 0
    metric::M = Chebyshev()
end

function estimate(infomeasure::CMI{Nothing}, e::Renyi, est::Rahimzamani, X, Y, Z)
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
        dmax = ds_joint[i]
        k̂ = dmax == 0 ? inrangecount(tree_joint, joint[i], 0.0) - 1  : k
        condmi += digamma(k̂)
        condmi -= log(inrangecount(tree_xz, XZ[i], dmax))
        condmi -= log(inrangecount(tree_yz, YZ[i], dmax))
        condmi += log(inrangecount(tree_z, Z[i], dmax))
    end
    condmi /= N

    return condmi / log(e.base, ℯ)
end
