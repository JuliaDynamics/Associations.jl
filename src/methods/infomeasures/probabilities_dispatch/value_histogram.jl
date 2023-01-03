# --------------------------------------------------------------------------------------
# Conditional mutual information
# --------------------------------------------------------------------------------------
function marginal_probabilities(
        measure::CMI{<:EntropyDefinition, <:CMIH4},
        est::ValueHistogram{<:FixedRectangularBinning{D}},
        x, y, z) where {D}
    e = measure.e
    X, Y, Z = Dataset(x), Dataset(y), Dataset(z);
    XZ = Dataset(X, Z)
    YZ = Dataset(Y, Z)
    XYZ = Dataset(X, Y, Z)
    DXZ, DYZ, DXYZ, DZ = dimension(XZ), dimension(YZ), dimension(XYZ), dimension(Z)

    # TODO: Figure out a clever, non-complicated way of using dimension-variable-width bins.
    # This is possible if encoding gives tuples of integers indicating the bins.
    # We can then just find the joint bin visitations, and then get the marginals from that.
    # But this approach is naively much slower, because we have to subset the joint many
    # timer, which gets expensive when input dimension and number of points is large.
    # Enforcing D == 1 ensures computes are reasonably quick.
    D == 1 || throw(ArgumentError("""Estimator ValueHistogram{<:FixedRectangularBinning{D}}\
     is only defined for `condmutualinfo` when `D == 1` (got D=$(D)). This is because our \
        implementation only handles hypersquare bins at the moment."""))
    ϵmin = est.binning.ϵmin[1]
    ϵmax = est.binning.ϵmax[1]
    N = est.binning.N
    estXZ = ValueHistogram(FixedRectangularBinning(ϵmin, ϵmax, N, DXZ))
    estYZ = ValueHistogram(FixedRectangularBinning(ϵmin, ϵmax, N, DYZ))
    estXYZ = ValueHistogram(FixedRectangularBinning(ϵmin, ϵmax, N, DXYZ))
    estZ = ValueHistogram(FixedRectangularBinning(ϵmin, ϵmax, N, DZ))

    pXZ = probabilities(estXZ, XZ)
    pYZ = probabilities(estYZ, YZ)
    pXYZ = probabilities(estXYZ, XYZ)
    pZ = probabilities(estZ, Z)
    return pXZ, pYZ, pXYZ, pZ
end



# Alternative approach for RectangularBinning, which ensures the same bins are used
# for all marginals, but is slower.
# function estimate(def::CMIH4, measure::CMI, est::VisitationFrequency{B}, x, y, z) where B
#     # Assumes `Dataset`-like structure of the input, i.e. rows are observations.
#     DX, DY, DZ = size(x, 2), size(y, 2), size(z, 2)
#     e = measure.e
#     b = est.binning

#     # We're going to encode the joint space, then subsample that to compute probabilities
#     # for the marginal entropies. Pre-compute these.
#     idxs_X = 1:DX
#     idxs_Y = (maximum(idxs_X) + 1):DX+DY
#     idxs_Z = (maximum(idxs_Y) + 1):DX+DY+DZ
#     idxs_XZ = [idxs_X; idxs_Z]
#     idxs_YZ = [idxs_Y; idxs_Z]

#     XYZ = Dataset(Dataset(x), Dataset(y), Dataset(z))
#     rb = RectangularBinEncoding(XYZ, b)
#     XYZ_bins = Dataset([decode(encode(rb, xyzᵢ)) for xyzᵢ in XYZ])
#     probs_xyz = probabilities(XYZ_bins)
#     probs_z = probabilities(XYZ_bins[:, idxs_Z])
#     probs_xz = probabilities(XYZ_bins[:, idxs_XZ])
#     probs_yz = probabilities(XYZ_bins[:, idxs_YZ])

#     condmutualinfo = entropy(e, probs_xz) +
#         entropy(e, probs_yz) -
#         entropy(e, probs_xyz) -
#         entropy(e, probs_z)
#     return condmutualinfo / log(e.base, ℯ)
# end

# --------------------------------------------------------------------------------------
# Mutual information
# --------------------------------------------------------------------------------------

# FixedRectangularBinning demands explicitly specifying the dimension, so need specialized
# dispatch.
function marginal_probabilities(
        measure::MutualInformation{<:EntropyDefinition, <:MIH3},
        est::ValueHistogram{<:FixedRectangularBinning{D}},
        x, y) where {D}
    e = measure.e
    X, Y = Dataset(x), Dataset(y); XY = Dataset(X, Y)
    DX, DY, DXY = dimension(X), dimension(Y), dimension(XY)
    D == 1 || throw(ArgumentError("""Estimator ValueHistogram{<:FixedRectangularBinning{D}}\
     must have dimension D == 1. Got D=$(D)."))"""))
    ϵmin = est.binning.ϵmin[1]
    ϵmax = est.binning.ϵmax[1]
    N = est.binning.N
    estx = ValueHistogram(FixedRectangularBinning(ϵmin, ϵmax, N, DX))
    esty = ValueHistogram(FixedRectangularBinning(ϵmin, ϵmax, N, DY))
    estxy = ValueHistogram(FixedRectangularBinning(ϵmin, ϵmax, N, DXY))
    pX = probabilities(estx, X)
    pY = probabilities(est, Y)
    pXY = probabilities(estxy, XY)
    return pX, pY, pXY
end
