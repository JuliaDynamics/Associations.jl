export CMI4H

"""
    CMI4HShannon <: ConditionalMutualInformationDefinition
    CMI4HShannon()

A directive to compute Shannon conditional mutual information as a sum of four entropy terms:

```math
I(X; Y | Z) = H^S(X, Z) + H^S(Y, z) - H^S(X, Y, Z) - H^S(Z),
```

where ``H^S(\\cdot)`` is the generalized entropy of the argument. The formula applies for
the continuous case too, just changing notation to lower-case ``h`` to indicate that the
*differential* entropy is used:

```math
I(X; Y | Z) = h^S(X, Z) + h^S(Y, z) - h^S(X, Y, Z) - h^G(Z).
```
"""
struct CMI4H <: ConditionalMutualInformationDefinition end

function estimate(def::CMI4H, measure::CMIShannon,
    est::Union{EntropyEstimator, ProbabilitiesEstimator}, x, y, z)
    e = measure.e
    X = Dataset(x)
    Y = Dataset(y)
    Z = Dataset(z)
    XZ = Dataset(X, Z)
    YZ = Dataset(Y, Z)
    XYZ = Dataset(X, Y, Z)
    c =  entropy(e, est, XZ) +
        entropy(e, est, YZ) -
        entropy(e, est, XYZ) -
        entropy(e, est, Z)
    return c / log(e.base, ℯ)
end

estimate(measure::CMIShannon, est::Union{ProbabilitiesEstimator, EntropyEstimator},
        x, y, z) = estimate(CMI4H(), measure, est, x, y, z)

###########################################################################################
# Specialized dispatch.
# ---------------------
# For some information measures and/or methods, specialized dispatch might be necessary
# or more efficient. Define these methods below.
###########################################################################################
function estimate(def::CMI4H, measure::CMIShannon,
        est::VisitationFrequency{B}, x, y, z) where B
    #@show "heyo"
    # Assumes `Dataset`-like structure of the input, i.e. rows are observations.
    DX, DY, DZ = size(x, 2), size(y, 2), size(z, 2)
    e = measure.e
    b = est.binning

    # We're going to encode the joint space, then subsample that to compute probabilities
    # for the marginal entropies. Pre-compute these.
    idxs_X = 1:DX
    idxs_Y = (maximum(idxs_X) + 1):DX+DY
    idxs_Z = (maximum(idxs_Y) + 1):DX+DY+DZ
    idxs_XZ = [idxs_X; idxs_Z]
    idxs_YZ = [idxs_Y; idxs_Z]

    XYZ = Dataset(Dataset(x), Dataset(y), Dataset(z))
    rb = RectangularBinEncoding(XYZ, b)
    XYZ_bins = Dataset([Entropies.encode_as_bin(xyzᵢ, rb) for xyzᵢ in XYZ])
    probs_xyz = probabilities(XYZ_bins)
    probs_z = probabilities(XYZ_bins[:, idxs_Z])
    probs_xz = probabilities(XYZ_bins[:, idxs_XZ])
    probs_yz = probabilities(XYZ_bins[:, idxs_YZ])

    cmi = entropy(e, probs_xz) +
        entropy(e, probs_yz) -
        entropy(e, probs_xyz) -
        entropy(e, probs_z)
    return cmi / log(e.base, ℯ)
end
