using ComplexityMeasures: SymbolicPermutation

# --------------------------------------------------------------------------------------
# Conditional mutual information
# --------------------------------------------------------------------------------------
function marginal_probabilities(
        measure::CMI,
        est::SymbolicPermutation{m},
        x::AbstractVector,
        y::AbstractVector,
        z::AbstractVector) where {m}
    πX, πY, πZ = marginal_encodings(est, x, y, z)
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
function marginal_probabilities(measure::MutualInformation,
        est::SymbolicPermutation{m}, x::AbstractVector, y::AbstractVector) where {m}
    πx, πy = marginal_encodings(est, x, y)
    πxy = Dataset(πx, πy)
    HX = probabilities(CountOccurrences(), πx)
    HY = probabilities(CountOccurrences(), πy)
    HXY = probabilities(CountOccurrences(), πxy)
    return HX, HY, HXY
end
