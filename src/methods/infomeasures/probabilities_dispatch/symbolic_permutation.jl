using ComplexityMeasures: SymbolicPermutation

# --------------------------------------------------------------------------------------
# Conditional mutual information
# --------------------------------------------------------------------------------------
function marginal_probabilities(
        measure::CMI{<:EntropyDefinition, <:CMIH4},
        est::SymbolicPermutation{m},
        x::AbstractVector,
        y::AbstractVector,
        z::AbstractVector) where {m}
    ex = embed(x, m, est.τ)
    ey = embed(y, m, est.τ)
    ez = embed(z, m, est.τ)

    πX = zeros(Int, length(ex))
    πY = zeros(Int, length(ey))
    πZ = zeros(Int, length(ez))
    @inbounds for (i, (x, y, z)) in enumerate(zip(ex, ey, ez))
        πX[i] = encode(est.encoding, x)
        πY[i] = encode(est.encoding, y)
        πZ[i] = encode(est.encoding, z)
    end
    πXZ = Dataset(πX, πZ)
    πYZ = Dataset(πY, πZ)
    πXYZ = Dataset(πX, πY, πZ)

    pXZ = probabilities(CountOccurrences(), πXZ)
    pYZ = probabilities(CountOccurrences(), πYZ)
    pXYZ = probabilities(CountOccurrences(), πXYZ)
    pZ = probabilities(CountOccurrences(), πZ)
    return pXZ, pYZ, pXYZ, pZ
end


# --------------------------------------------------------------------------------------
# Mutual information
# --------------------------------------------------------------------------------------
 # TODO: this can be extended to any symbolic weighted scheme, but then we need to
# explictly incorporate the weights in the encodings. I have a method paper in prep about
# the AmplitudeAwareTransferEntropy, which this fits nicely in with.
function marginal_probabilities(measure::MutualInformation{<:EntropyDefinition, <:MIH3},
        est::SymbolicPermutation{m}, x::AbstractVector, y::AbstractVector) where {m}
    e = measure.e
    ex = embed(x, m, est.τ)
    ey = embed(y, m, est.τ)
    πx = zeros(Int, length(ex))
    πy = zeros(Int, length(ey))
    @inbounds for (i, (x, y)) in enumerate(zip(ex, ey))
        πx[i] = encode(est.encoding, x)
        πy[i] = encode(est.encoding, y)
    end
    πxy = Dataset(πx, πy)

    HX = probabilities(CountOccurrences(), πx)
    HY = probabilities(CountOccurrences(), πy)
    HXY = probabilities(CountOccurrences(), πxy)
    return HX, HY, HXY
end
