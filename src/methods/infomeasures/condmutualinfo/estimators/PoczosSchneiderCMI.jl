using ComplexityMeasures: ball_volume
using SpecialFunctions: gamma

export PoczosSchneiderCMI
"""
    PoczosSchneiderCMI <: ConditionalMutualInformationEstimator
    PoczosSchneiderCMI(k = 1, w = 0)

The `PoczosSchneiderCMI` estimator computes various (differential) conditional
mutual informations, using a `k`-th nearest neighbor approach (Póczos & Schneider,
2012)[^Póczos2012].

## Description

Póczos & Schneider (2012) defines the following quantities.

### Rényi conditional mutual information

[^Póczos2012]:
    Póczos, B., & Schneider, J. (2012, March). Nonparametric estimation of conditional
    information and divergences. In Artificial Intelligence and Statistics (pp. 914-923).
    PMLR.
"""
Base.@kwdef struct PoczosSchneiderCMI{M} <: ConditionalMutualInformationEstimator
    k::Int = 1
    w::Int = 0
    metric::M = Euclidean()
end

function estimate(measure::CMIRenyi, est::PoczosSchneiderCMI, x, y, z)
    e = measure.e
    c = 1/(e.q-1)*log(Q3(e, est, x, y, z))
    return c / log(e.base, ℯ)
end

function Q3(e::EntropyDefinition, est::PoczosSchneiderCMI, x, y, z)
    q = e.q
    (; k, w, metric) = est
    @assert length(x) == length(y) == length(z)
    N = length(x)
    YZ = Dataset(y, z)
    XZ = Dataset(x, z)
    XYZ = Dataset(x, y, z)
    Z = Dataset(z)

    tree_YZ = KDTree(YZ, metric)
    tree_XZ = KDTree(XZ, metric)
    tree_XYZ = KDTree(XYZ, metric)
    tree_Z = KDTree(Z, metric)

    idxs_YZ, dists_YZ = bulksearch(tree_YZ, YZ, NeighborNumber(k), Theiler(w))
    idxs_XZ, dists_XZ = bulksearch(tree_XZ, XZ, NeighborNumber(k), Theiler(w))
    idxs_XYZ, dists_XYZ = bulksearch(tree_XYZ, XYZ, NeighborNumber(k), Theiler(w))
    idxs_Z, dists_Z = bulksearch(tree_Z, Z, NeighborNumber(k), Theiler(w))

    ds_YZ = last.(dists_YZ)
    ds_XZ = last.(dists_XZ)
    ds_XYZ = last.(dists_XYZ)
    ds_z = last.(dists_Z)

    # Not sure about the index sets here.
    # fyz = (N - 1)^(1 - q)
    #fxz = (N - 1)^(1 - q)
    bv_yz = ball_volume(dimension(YZ))
    bv_xz = ball_volume(dimension(XZ))
    bv_xyz = ball_volume(dimension(XYZ))
    bv_z = ball_volume(dimension(Z))

    # The ninality of each set is the same, because all variables are equal-length
    n = (N - 1)
    f = ( (bv_xyz * n)^(1 - q) * (bv_z  * n)^(1 - q) ) /
        ( (bv_xz  * n)^(1 - q) * (bv_yz * n)^(1 - q) )
    B² = (gamma(k)^2 / (gamma(k - q + 1)*gamma(k + q - 1)))^2
    condmi = 0.0
    for i = 1:N
        condmi += f * (ds_XYZ[i] * ds_z[i]) / (ds_XZ[i] * ds_YZ[i]) * B²
    end
    condmi /= N
    return condmi
end
