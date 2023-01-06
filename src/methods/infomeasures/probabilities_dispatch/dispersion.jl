
# --------------------------------------------------------------------------------------
# Conditional mutual information
# --------------------------------------------------------------------------------------

# This will work for any future versions of `Dispersion` with different encoding schemes.
function marginal_probabilities(
        measure::CMI,
        est::Dispersion,
        x::AbstractVector, y::AbstractVector, z::AbstractVector)
    e = measure.e
    # These are now symbol time series
    πX, πY, πZ = marginal_encodings(est, x, y, z)
    πXZ = Dataset(πX, πZ)
    πYZ = Dataset(πY, πZ)
    πXYZ = Dataset(πX, πY, πZ)

    if est.m == 1 # d == 1 is the same as not performing any embedding.
        pXZ = probabilities(CountOccurrences(), πXZ)
        pYZ = probabilities(CountOccurrences(), πYZ)
        pXYZ = probabilities(CountOccurrences(), πXYZ)
        pZ = probabilities(CountOccurrences(), πZ)
        return pXZ, pYZ, pXYZ, pZ
    else
        # To make this analogous to the dispersion entropy, for `d > 1`, embed each
        # symbol timeseries and compute histograms over the embeddings.
        m, τ = est.m, est.τ
        τs = tuple((x for x in 0:-τ:-(m-1)*τ)...)
        EπX = genembed(πX, τs, ones(m))
        EπY = genembed(πY, τs, ones(m))
        EπZ = genembed(πZ, τs, ones(m))
        EπXZ = Dataset(EπX, EπZ)
        EπYZ = Dataset(EπY, EπZ)
        EπXYZ = Dataset(EπX, EπY, EπZ)

        pXZ = probabilities(CountOccurrences(), EπXZ)
        pYZ = probabilities(CountOccurrences(), EπYZ)
        pXYZ = probabilities(CountOccurrences(), EπXYZ)
        pZ = probabilities(CountOccurrences(), EπZ)
        return pXZ, pYZ, pXYZ, pZ
    end
end


# --------------------------------------------------------------------------------------
# Mutual information
# --------------------------------------------------------------------------------------

# This will work for any future versions of `Dispersion` with different encoding schemes.
function marginal_probabilities(measure::MutualInformation,
        est::Dispersion, x::AbstractVector, y::AbstractVector)
    e = measure.e
    # These are now symbol time series
    πx, πy = marginal_encodings(est, x, y)
    πxy = Dataset(πx, πy)

    if est.m == 1 # d == 1 is the same as not performing any embedding.
        pX = probabilities(CountOccurrences(), πx)
        pY = probabilities(CountOccurrences(), πy)
        pXY = probabilities(CountOccurrences(), πxy)
        return pX, pY, pXY
    else
        # To make this analogous to the dispersion entropy, for `d > 1`, embed each
        # symbol timeseries and compute histograms over the embeddings.
        m, τ = est.m, est.τ
        τs = tuple((x for x in 0:-τ:-(m-1)*τ)...)
        Eπx = genembed(πx, τs, ones(m))
        Eπy = genembed(πy, τs, ones(m))
        Eπxy = Dataset(Eπx, Eπy)
        pX = probabilities(CountOccurrences(), Eπx)
        pY = probabilities(CountOccurrences(), Eπy)
        pXY = probabilities(CountOccurrences(), Eπxy)
        return pX, pY, pXY
    end
end
