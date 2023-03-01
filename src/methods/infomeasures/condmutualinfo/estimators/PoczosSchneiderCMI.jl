using ComplexityMeasures: ball_volume
using SpecialFunctions: gamma

export PoczosSchneiderCMI
"""
    PoczosSchneiderCMI <: ConditionalMutualInformationEstimator
    PoczosSchneiderCMI(k = 1, w = 0)

The `PoczosSchneiderCMI` estimator computes various (differential) conditional
mutual informations, using a `k`-th nearest neighbor approach (Póczos & Schneider,
2012)[^Póczos2012].

[^Póczos2012]:
    Póczos, B., & Schneider, J. (2012, March). Nonparametric estimation of conditional
    information and divergences. In Artificial Intelligence and Statistics (pp. 914-923).
    PMLR.
"""
Base.@kwdef struct PoczosSchneiderCMI{M} <: ConditionalMutualInformationEstimator
    k::Int = 1
    w::Int = 0
    metric::M = Euclidean() # Needs to be euclidean for ball volume formula to be valid.
end

function estimate(measure::CMIRenyiPoczos, est::PoczosSchneiderCMI, x, y, z)
    e = measure.e
    c = log(Q3(e, est, x, y, z)) / (e.q-1)
    return c / log(ℯ, e.base)
end

function Q3(e::EntropyDefinition, est::PoczosSchneiderCMI, x, y, z)
    q = e.q
    (; k, w, metric) = est
    @assert length(x) == length(y) == length(z)
    N = length(x)
    YZ = StateSpaceSet(y, z)
    XZ = StateSpaceSet(x, z)
    XYZ = StateSpaceSet(x, y, z)
    Z = StateSpaceSet(z)

    tree_YZ = KDTree(YZ, metric)
    tree_XZ = KDTree(XZ, metric)
    tree_XYZ = KDTree(XYZ, metric)
    tree_Z = KDTree(Z, metric)

    idxs_YZ, dists_YZ = bulksearch(tree_YZ, YZ, NeighborNumber(k), Theiler(w))
    idxs_XZ, dists_XZ = bulksearch(tree_XZ, XZ, NeighborNumber(k), Theiler(w))
    idxs_XYZ, dists_XYZ = bulksearch(tree_XYZ, XYZ, NeighborNumber(k), Theiler(w))
    idxs_Z, dists_Z = bulksearch(tree_Z, Z, NeighborNumber(k), Theiler(w))

    ds_YZ = last.(dists_YZ) .^ (dimension(YZ) * (1 - q))
    ds_XZ = last.(dists_XZ) .^ (dimension(XZ) * (1 - q))
    ds_XYZ = last.(dists_XYZ) .^ (dimension(XYZ) * (1 - q))
    ds_z = last.(dists_Z) .^ (dimension(Z) * (1 - q))

    # Not sure about the index sets here.
    # fyz = (N - 1)^(1 - q)
    #fxz = (N - 1)^(1 - q)
    bv_yz = ball_volume(dimension(YZ)) ^ (1 - q)
    bv_xz = ball_volume(dimension(XZ)) ^ (1 - q)
    bv_xyz = ball_volume(dimension(XYZ)) ^(1 - q)
    bv_z = ball_volume(dimension(Z)) ^ (1 - q)

    # The cardinality of each set is the same, because all variables are equal-length
    n = (N - 1)
    f = ( (bv_xyz * n) * (bv_z  * n) ) / ( (bv_xz  * n) * (bv_yz * n) )
    B² = (gamma(k)^2 / (gamma(k - q + 1)*gamma(k + q - 1)))^2
    condmi = 0.0
    for i = 1:N
        condmi += f * (ds_XYZ[i] * ds_z[i]) / (ds_XZ[i] * ds_YZ[i]) * B²
    end
    condmi /= N
    return condmi
end


function condmutualinfo(est::PoczosSchneiderCMI, x, y, z)
    return condmutualinfo(CMIRenyiPoczos(), est, x, y, z)
end
